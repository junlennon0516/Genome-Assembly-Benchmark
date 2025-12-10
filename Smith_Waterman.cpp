#include <iostream>
#include <vector>
#include <string>
#include <string_view>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <iomanip>

struct ScoringSystem {
    int matchScore = 3;
    int mismatchPenalty = -2;
    int gapPenalty = -3;
};

class SmithWatermanAligner {
public:
    SmithWatermanAligner(ScoringSystem scores) : scores_(scores) {}

    // 반환값: 점수, 최대 점수 위치(참조 게놈 상의 인덱스)
    std::pair<int, int> align(std::string_view reference, std::string_view read, int maxMismatch) {
        size_t refLen = reference.size();
        size_t readLen = read.size();

        // 메모리 재할당 방지를 위해 클래스 멤버 혹은 static 버퍼 사용 고려 가능하나
        // 여기서는 vector의 오버헤드를 줄이기 위해 지역 변수 재사용
        std::vector<int> prevRow(readLen + 1, 0);
        std::vector<int> currRow(readLen + 1, 0);

        int maxScore = 0;
        int maxPos = -1;

        for (size_t i = 1; i <= refLen; ++i) {
            int mismatchCount = 0; 
            
            for (size_t j = 1; j <= readLen; ++j) {
                // 1. Match/Mismatch 계산
                int matchCalc;
                bool isMatch = reference[i - 1] == read[j - 1];
                
                if (!isMatch) mismatchCount++; // 단순화된 mismatch 카운팅 (엄밀한 경로 추적은 백트래킹 필요)
                
                // Mismatch 한계 초과 시 패널티 부여 방식 조정
                if (!isMatch && mismatchCount > maxMismatch) {
                   matchCalc = -9999; // 강제 경로 차단
                } else {
                   matchCalc = prevRow[j - 1] + (isMatch ? scores_.matchScore : scores_.mismatchPenalty);
                }

                // 2. Gap 계산
                int del = prevRow[j] + scores_.gapPenalty;
                int ins = currRow[j - 1] + scores_.gapPenalty;

                // 3. Max Score (Local Alignment는 음수일 경우 0)
                currRow[j] = std::max({0, matchCalc, del, ins});

                if (currRow[j] > maxScore) {
                    maxScore = currRow[j];
                    maxPos = i - readLen; // 정렬 시작 위치 추정
                }
            }
            // Swap rows
            prevRow = currRow; // Vector copy 대신 move/swap이 효율적이나 vector<int>는 작으므로 swap
            std::fill(currRow.begin(), currRow.end(), 0);
        }
        return {maxScore, maxPos};
    }

private:
    ScoringSystem scores_;
};

class GenomeProcessor {
public:
    static void run(const std::string& refFile, const std::string& readsFile, const std::string& outFile) {
        std::string reference = readFile(refFile);
        std::string myGenome = readFile("My_genome.txt"); // 검증용
        std::vector<std::string> reads;
        
        {
            std::ifstream file(readsFile);
            std::string line;
            while (file >> line) reads.push_back(line);
        }

        SmithWatermanAligner aligner({3, -2, -3});
        std::vector<std::pair<int, std::string>> restoredReads;
        int maxMismatch = 2;

        std::cout << "Processing " << reads.size() << " reads..." << std::endl;
        auto start = std::chrono::high_resolution_clock::now();

        for (const auto& read : reads) {
            auto result = aligner.align(reference, read, maxMismatch);
            int score = result.first;
            int pos = result.second;

            if (pos >= 0 && pos + read.size() <= reference.size() && score > 0) {
                restoredReads.emplace_back(pos, read);
            }
        }

        // 위치 기반 정렬
        std::sort(restoredReads.begin(), restoredReads.end());

        // 복원 (Reconstruction)
        std::string restoredGenome(reference.size(), 'N');
        for (const auto& p : restoredReads) {
            int startPos = p.first;
            const std::string& seq = p.second;
            for (size_t i = 0; i < seq.size(); ++i) {
                if (startPos + i < restoredGenome.size() && restoredGenome[startPos + i] == 'N') {
                    restoredGenome[startPos + i] = seq[i];
                }
            }
        }

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "Sequencing took: " << elapsed.count() << " seconds" << std::endl;

        std::ofstream out(outFile);
        out << restoredGenome;
        
        verifyResults(myGenome, restoredGenome);
    }

private:
    static std::string readFile(const std::string& filename) {
        std::ifstream file(filename);
        if (!file) { std::cerr << "Open failed: " << filename << std::endl; exit(1); }
        std::string str; 
        std::getline(file, str); 
        return str;
    }

    static void verifyResults(std::string_view original, std::string_view assembled) {
        size_t matches = 0;
        size_t len = std::min(original.size(), assembled.size());
        for (size_t i = 0; i < len; ++i) {
            if (original[i] == assembled[i]) matches++;
        }
        std::cout << "Accuracy: " << (double)matches / original.size() * 100.0 << "%" << std::endl;
    }
};

int main() {
    GenomeProcessor::run("reference.txt", "reads.txt", "restored.txt");
    return 0;
}