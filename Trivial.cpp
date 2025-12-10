#include <iostream>
#include <vector>
#include <string>
#include <string_view>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <map>
#include <iomanip>

class ReferenceAssembler {
public:
    ReferenceAssembler(const std::string& refFile, const std::string& readsFile) 
        : maxMismatches_(2) {
        loadGenome(refFile);
        loadReads(readsFile);
    }

    void runAssembly(const std::string& outputFile) {
        std::string assembledGenome(referenceGenome_.size(), 'N');
        // 각 위치별 염기 빈도수 저장: [position][A,C,G,T,N] -> map 대신 vector나 array가 빠름
        // 하지만 편의상 map<int, map<char, int>> 구조 유지하되 position 접근 최적화
        std::map<int, std::map<char, int>> baseCounts; 

        std::cout << "Mapping reads..." << std::endl;
        auto start = std::chrono::high_resolution_clock::now();

        size_t totalReads = reads_.size();
        for (size_t i = 0; i < totalReads; ++i) {
            auto match = findBestMatch(reads_[i]);
            if (match.position != -1) {
                recordRead(match, reads_[i], baseCounts);
            }
            if ((i + 1) % (totalReads / 100 + 1) == 0) printProgress(i + 1, totalReads);
        }

        buildConsensus(assembledGenome, baseCounts);
        
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "\nTime: " << std::chrono::duration<double>(end - start).count() << "s" << std::endl;

        saveFile(outputFile, assembledGenome);
        verify(assembledGenome);
    }

private:
    std::string referenceGenome_;
    std::vector<std::string> reads_;
    const int maxMismatches_;

    struct MatchResult {
        int position = -1;
        std::vector<int> mismatchIndices; // read 내부 인덱스
    };

    void loadGenome(const std::string& f) {
        std::ifstream file(f);
        if(!file) throw std::runtime_error("No ref file");
        std::getline(file, referenceGenome_);
    }

    void loadReads(const std::string& f) {
        std::ifstream file(f);
        if(!file) throw std::runtime_error("No reads file");
        std::string s;
        while(std::getline(file, s)) reads_.push_back(s);
    }

    MatchResult findBestMatch(std::string_view read) {
        MatchResult best;
        int minMismatches = maxMismatches_ + 1;
        int readLen = read.length();
        int maxIdx = referenceGenome_.length() - readLen;

        for (int i = 0; i <= maxIdx; ++i) {
            int currentMismatches = 0;
            std::vector<int> mismatches;
            
            // 성능 핵심 구간 (Naive 비교)
            bool possible = true;
            for (int j = 0; j < readLen; ++j) {
                if (referenceGenome_[i + j] != read[j]) {
                    currentMismatches++;
                    if (currentMismatches > maxMismatches_) {
                        possible = false;
                        break;
                    }
                    mismatches.push_back(j);
                }
            }

            if (possible && currentMismatches < minMismatches) {
                minMismatches = currentMismatches;
                best.position = i;
                best.mismatchIndices = std::move(mismatches);
                if (minMismatches == 0) break; // 완전 일치 찾으면 조기 종료
            }
        }
        return best;
    }

    void recordRead(const MatchResult& match, std::string_view read, std::map<int, std::map<char, int>>& counts) {
        int startPos = match.position;
        // 불일치 위치 Set으로 변환하여 빠른 조회
        std::vector<bool> isMismatch(read.size(), false);
        for(int idx : match.mismatchIndices) isMismatch[idx] = true;

        for(size_t i=0; i<read.size(); ++i) {
            if (startPos + i < referenceGenome_.size()) {
                // 원본 소스의 로직: "불일치 위치가 아니면 Read의 염기를, 불일치면...?"
                // 원본 로직 유지: Mismatch 리스트에 있는 위치는 기록에 포함시킴(read의 염기로)
                // 하지만 보통 Assembly에서는 read 염기를 신뢰하여 투표함.
                counts[startPos + i][read[i]]++;
            }
        }
    }

    void buildConsensus(std::string& assembled, const std::map<int, std::map<char, int>>& counts) {
        for (auto const& [pos, bases] : counts) {
            char bestBase = 'N';
            int maxCount = -1;
            for (auto const& [base, count] : bases) {
                if (count > maxCount) {
                    maxCount = count;
                    bestBase = base;
                }
            }
            if (pos < assembled.size()) assembled[pos] = bestBase;
        }
    }

    void printProgress(size_t done, size_t total) {
        double p = (double)done / total * 100.0;
        std::cout << "\rProgress: " << std::fixed << std::setprecision(1) << p << "%" << std::flush;
    }

    void saveFile(const std::string& f, const std::string& data) {
        std::ofstream(f) << data;
    }

    void verify(std::string_view assembled) {
        std::ifstream file("my_genome.txt");
        std::string original;
        std::getline(file, original);
        
        int matches = 0;
        size_t len = std::min(original.size(), assembled.size());
        for(size_t i=0; i<len; ++i) if(original[i] == assembled[i]) matches++;
        
        std::cout << "Accuracy: " << (double)matches / original.size() * 100.0 << "%\n";
    }
};

int main() {
    try {
        ReferenceAssembler assembler("reference.txt", "reads.txt");
        assembler.runAssembly("assembled_genome.txt");
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
    return 0;
}