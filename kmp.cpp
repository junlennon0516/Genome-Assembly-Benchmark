#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <string_view>
#include <algorithm>
#include <unordered_set>
#include <cmath>
#include <chrono>
#include <iomanip>

// KMP 알고리즘을 담당하는 클래스
class KMPAligner {
public:
    static int match(std::string_view text, std::string_view pattern, int maxMismatch) {
        if (pattern.size() > text.size()) return -1;

        std::vector<int> lps = computeLPS(pattern);
        int i = 0, j = 0;
        int mismatchCount = 0;
        int n = text.size();
        int m = pattern.size();

        while (i < n) {
            if (j < m && text[i] == pattern[j]) {
                i++; j++;
            } else if (j < m && text[i] != pattern[j]) {
                mismatchCount++;
                if (mismatchCount > maxMismatch) {
                    mismatchCount = 0;
                    if (j != 0) j = lps[j - 1];
                    else i++;
                } else {
                    i++; j++;
                }
            }

            if (j == m && mismatchCount <= maxMismatch) {
                return i - j;
            }
        }
        return -1;
    }

private:
    static std::vector<int> computeLPS(std::string_view pattern) {
        int m = pattern.size();
        std::vector<int> lps(m, 0);
        int j = 0;
        for (int i = 1; i < m; i++) {
            while (j > 0 && pattern[i] != pattern[j]) {
                j = lps[j - 1];
            }
            if (pattern[i] == pattern[j]) {
                j++;
            }
            lps[i] = j;
        }
        return lps;
    }
};

// 파일 유틸리티
class FileUtils {
public:
    static std::string readFile(const std::string& filename) {
        std::ifstream file(filename);
        if (!file) throw std::runtime_error("Cannot open file: " + filename);
        return std::string((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    }

    static std::vector<std::string> readReads(const std::string& filename) {
        std::ifstream file(filename);
        if (!file) throw std::runtime_error("Cannot open file: " + filename);
        std::vector<std::string> reads;
        std::string read;
        while (file >> read) reads.push_back(read);
        return reads;
    }

    static void writeFile(const std::string& filename, const std::string& data) {
        std::ofstream file(filename);
        if (!file) throw std::runtime_error("Cannot write file: " + filename);
        file << data;
    }
};

// 메인 로직 처리
class GenomeAssembler {
public:
    static std::string mergeReads(const std::string& reference, const std::vector<std::string>& reads, int maxMismatch) {
        std::string merged = reference;
        std::unordered_set<std::string> processedReads;
        size_t totalReads = reads.size();
        
        std::cout << "Reads 병합 진행 중..." << std::endl;
        auto start = std::chrono::steady_clock::now();

        for (size_t i = 0; i < totalReads; ++i) {
            const auto& read = reads[i];

            // 이미 처리된 Read는 건너뜀 (Set 검색 최적화)
            if (processedReads.count(read)) continue;

            int pos = KMPAligner::match(merged, read, maxMismatch);
            if (pos != -1) {
                applyMerge(merged, read, pos);
                processedReads.insert(read);
            }

            // 진행률 표시 (1% 단위)
            if (totalReads > 100 && (i + 1) % (totalReads / 100) == 0) {
                std::cout << "\rProgress: " << ((i + 1) * 100) / totalReads << "% completed" << std::flush;
            }
        }

        auto end = std::chrono::steady_clock::now();
        std::cout << "\nTotal execution time: " 
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0 
                  << " sec\n";

        return merged;
    }

    static double calculateJaccardSimilarity(const std::string& seq1, const std::string& seq2, int k) {
        if (seq1.size() < k || seq2.size() < k) return 0.0;

        std::unordered_set<std::string> set1;
        for (size_t i = 0; i <= seq1.size() - k; ++i) set1.insert(seq1.substr(i, k));

        std::unordered_set<std::string> set2;
        for (size_t i = 0; i <= seq2.size() - k; ++i) set2.insert(seq2.substr(i, k));

        size_t intersection = 0;
        for (const auto& s : set1) {
            if (set2.count(s)) intersection++;
        }

        size_t unionSize = set1.size() + set2.size() - intersection;
        return unionSize == 0 ? 0.0 : static_cast<double>(intersection) / unionSize;
    }

private:
    static void applyMerge(std::string& target, std::string_view read, int pos) {
        for (size_t j = 0; j < read.size(); ++j) {
            size_t idx = pos + j;
            if (idx < target.size()) {
                if (target[idx] != read[j]) target[idx] = std::tolower(read[j]); // 불일치 마킹
            } else {
                target += std::tolower(read[j]);
            }
        }
    }
};

int main() {
    try {
        const std::string REF_FILE = "reference.txt";
        const std::string READS_FILE = "reads.txt";
        const std::string OUTPUT_FILE = "merged_genome.txt";
        const std::string MY_GENOME_FILE = "My_genome.txt";

        std::string reference = FileUtils::readFile(REF_FILE);
        std::vector<std::string> reads = FileUtils::readReads(READS_FILE);
        std::string myGenome = FileUtils::readFile(MY_GENOME_FILE);

        std::string mergedGenome = GenomeAssembler::mergeReads(reference, reads, 2);
        FileUtils::writeFile(OUTPUT_FILE, mergedGenome);

        int k = std::max(1, static_cast<int>(std::log2(myGenome.size()) / 2));
        double similarity = GenomeAssembler::calculateJaccardSimilarity(myGenome, mergedGenome, k);

        std::cout << "Jaccard Similarity (k=" << k << "): " << similarity * 100 << "%" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}