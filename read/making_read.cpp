#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;

// Reads를 스트림으로 생성하여 저장하는 함수
void generate_Reads(const string& input_file, const string& output_file, size_t read_length, size_t overlap)
{
    ifstream input(input_file);
    ofstream output(output_file);

    if (!input)
    {
        cerr << "파일을 열 수 없습니다." << input_file << endl;
        return;
    }
    if (!output)
    {
        cerr << "파일을 열 수 없습니다." << output_file << endl;
        return;
    }

    size_t step = read_length - overlap; // Reads 간 이동 거리
    string buffer; // 데이터를 임시로 저장할 버퍼
    buffer.reserve(read_length); // 버퍼 크기 예약

    size_t current_position = 0;
    size_t reads_count = 0; // Reads 개수 카운트
    char base;

    while (input.get(base))
    {
        // 버퍼에 현재 읽은 염기를 추가
        buffer += base;
        current_position++;

        // 버퍼가 Read 길이만큼 채워졌으면 Reads 생성
        if (buffer.size() == read_length)
        {
            output << buffer << endl; // Read 저장
            // 버퍼를 겹치는 부분만 남기고 초기화
            buffer = buffer.substr(step);
            reads_count++;
        }
    }

    input.close();
    output.close();

    cout << "Reads가 " << output_file << " 파일에 저장되었습니다." << endl;
}

int main()
{
    string input_file = "My_genome.txt";  // 변형된 유전체 파일
    string output_file = "reads.txt";    // 생성된 Reads 저장 파일

    // coverage 계산하는것도 필요
    size_t read_length = 60; // 각 Read의 길이
    size_t overlap = 56;      // Reads 간 겹치는 길이

    generate_Reads(input_file, output_file, read_length, overlap);

    return 0;
}