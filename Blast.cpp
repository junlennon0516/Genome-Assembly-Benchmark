#include <iostream>
#include <vector>
#include <string>
#include <string_view>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <map>
#include <unordered_map>
#include <iomanip>
#include <limits>

class BlastAssembler {
public:
    BlastAssembler(const std::string& refFile, const std::string& readsFile) 
        : wordSize_(11), maxMismatches_(2) {
        loadData(refFile, readsFile);
        createIndex();
    }

    void run(const std::string& outFile) {
        std::string assembled(reference_.size(), 'N');
        std::map<int, std::map<char, int>> consensusMap;

        std::cout << "BLAST Assembly Start..." << std::endl;
        auto start = std::chrono::high_resolution_clock::now();

        size_t total = reads_.size();
        for(size_t i=0; i<total; ++i) {
            auto result = findBestMatch(reads_[i]);
            if (result.first != -1) {
                updateConsensus(consensusMap, result.first, reads_[i]);
            }
            if ((i+1) % (total/100 + 1) == 0) 
                std::cout << "\rProgress: " << ((i+1)*100)/total << "%" << std::flush;
        }

        finalizeGenome(assembled, consensusMap);
        
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "\nTime: " << std::chrono::duration<double>(end - start).count() << "s\n";
        
        saveAndVerify(outFile, assembled);
    }

private:
    std::string reference_;
    std::vector<std::string> reads_;
    std::unordered_map<std::string, std::vector<int>> kmerIndex_;
    const int wordSize_;
    const int maxMismatches_;

    void loadData(const std::string& rFile, const std::string& qFile) {
        std::ifstream rf(rFile), qf(qFile);
        if(!rf || !qf) throw std::runtime_error("File error");
        std::getline(rf, reference_);
        std::string line;
        while(std::getline(qf, line)) reads_.push_back(line);
    }

    void createIndex() {
        if (reference_.size() < wordSize_) return;
        for (size_t i = 0; i <= reference_.size() - wordSize_; ++i) {
            kmerIndex_[reference_.substr(i, wordSize_)].push_back(i);
        }
    }

    std::pair<int, int> findBestMatch(const std::string& query) {
        // 1. Find HSP Candidates
        std::map<int, int> hits; // offset -> count
        for (size_t i = 0; i <= query.size() - wordSize_; ++i) {
            std::string word = query.substr(i, wordSize_);
            if (kmerIndex_.count(word)) {
                for (int refPos : kmerIndex_[word]) {
                    int startPos = refPos - i;
                    if (startPos >= 0) hits[startPos]++;
                }
            }
        }

        // 2. Extend & Score
        int bestPos = -1;
        int maxScore = std::numeric_limits<int>::min();

        for (auto const& [pos, count] : hits) {
            if (count < 1) continue; // 최소 1개 이상의 word match 필요 (원 코드는 2였으나 짧은 read 고려)
            
            // Extension (검증)
            int score = 0;
            int mismatches = 0;
            bool valid = true;
            
            for (size_t j = 0; j < query.size(); ++j) {
                if (pos + j >= reference_.size()) { valid = false; break; }
                
                if (reference_[pos + j] == query[j]) {
                    score += 2;
                } else {
                    score -= 1;
                    mismatches++;
                }
            }

            if (valid && mismatches <= maxMismatches_) {
                if (score > maxScore) {
                    maxScore = score;
                    bestPos = pos;
                }
            }
        }
        return {bestPos, maxScore};
    }

    void updateConsensus(std::map<int, std::map<char, int>>& map, int pos, const std::string& read) {
        for (size_t i = 0; i < read.size(); ++i) {
            map[pos + i][read[i]]++;
        }
    }

    void finalizeGenome(std::string& genome, const std::map<int, std::map<char, int>>& map) {
        for (auto const& [idx, counts] : map) {
            if (idx >= genome.size()) continue;
            char best = 'N';
            int maxC = -1;
            for (auto const& [base, cnt] : counts) {
                if (cnt > maxC) { maxC = cnt; best = base; }
            }
            genome[idx] = best;
        }
    }

    void saveAndVerify(const std::string& filename, const std::string& assembled) {
        std::ofstream(filename) << assembled;
        
        std::ifstream f("My_genome.txt");
        if(f) {
            std::string original; std::getline(f, original);
            int correct = 0;
            size_t len = std::min(original.size(), assembled.size());
            for(size_t i=0; i<len; ++i) if(original[i] == assembled[i]) correct++;
            std::cout << "Accuracy: " << (double)correct/original.size()*100.0 << "%\n";
        }
    }
};

int main() {
    try {
        BlastAssembler ba("reference.txt", "reads.txt");
        ba.run("assembled_genome.txt");
    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
    }
    return 0;
}