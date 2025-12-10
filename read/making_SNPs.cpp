#include <iostream>
#include <fstream>
#include <vector>
#include <random>
using namespace std;

// SNP 위치 생성 함수(1000bp당 1개, 위치는 1~1000사이 랜덤)
vector<size_t> generate_SNP_position(size_t length)
{
    vector<size_t> SNP_position;
    size_t SNP_total = length * 0.05; // SNP 생성 비율 설정
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> dist(20.0, 4.0); // 평균 20bp, 표준편차 4bp 
    // -> SNP 간격이 생물학적으로 더 현실적인 분포를 가지도록 조정

    size_t current_position = 0;
    for (size_t i = 0; i < SNP_total; ++i)
    {
        // 다음 SNP까지의 거리 생성 (최소 1bp 이상)
        int distance_to_next_SNP = max(1, static_cast<int>(dist(gen)));
        current_position += distance_to_next_SNP;

        // 시퀀스를 초과하지 않으면 SNP 위치 저장
        if (current_position < length)
        {
            SNP_position.push_back(current_position);
        }
        else
        {
            break; // 시퀀스 끝에 도달하면 중단
        }
    }

    return SNP_position;
}

// 염기를 랜덤하게 변형하는 함수
char transform_base(char base)
{
    string bases = "AGCT";
    bases.erase(bases.find(base), 1); // 원래 염기를 제외
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dist(0, 2);

    return bases[dist(gen)];
}

int main()
{
    string reference_file = "reference.txt"; // 원본 파일
    string output_file = "My_genome.txt";    // 변형된 염기서열 저장 파일

    // SNP 위치 생성
    const size_t length = 10000; // 예상 염기서열 길이
    vector<size_t> SNP_position = generate_SNP_position(length);

    // 파일 스트림 열기
    ifstream input(reference_file);
    ofstream output(output_file);

    if (!input)
    {
        cerr << "파일을 열 수 없습니다." << reference_file << endl;
        return 1;
    }
    if (!output)
    {
        cerr << "파일을 열 수 없습니다." << output_file << endl;
        return 1;
    }

    // 염기서열 처리 (스트림 방식)
    char base;
    size_t position = 0;         // 현재 파일에서 읽은 위치
    size_t SNP_index = 0;        // SNP 위치 벡터의 현재 인덱스

    // 한 글자씩 읽기
    while (input.get(base))
    {
        if (SNP_index < SNP_position.size() && position == SNP_position[SNP_index])
        {
            char original_base = base; //-> 바뀐 위치 출력해보는 코드
            base = transform_base(base);    // SNP 위치에서 염기 변형
            ++SNP_index;                    // 다음 SNP로 이동
        }
        output.put(base);                   // 변형된 염기를 출력 파일에 저장
        ++position;                         // 현재 위치 업데이트
    }

    input.close();
    output.close();

    cout << "변형된 SNP 위치의 개수: " << SNP_position.size() << endl;
    cout << "변형된 염기서열이 " << output_file << " 파일에 저장되었습니다." << endl;

    return 0;
}