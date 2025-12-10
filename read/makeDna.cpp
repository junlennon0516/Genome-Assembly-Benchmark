#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <ctime>

#define MAX_NUM 10000 // 생성할 DNA 시퀀스의 총 길이

using namespace std;
int main() {
	string line;
	ofstream file("reference.txt"); // 출력 파일 스트림 생성
	// DNA 염기 서열을 표현하는 4가지 문자
	char nucleic_sequnce[4] = { 'A', 'C', 'G', 'T' };
	int random_num;

	random_device rd; // random_device를 사용하여 시드 생성
	mt19937 gen(rd()); // 난수 생성기 초기화
	//random number의 범위: 0에서부터 3까지
	uniform_int_distribution<int> dis(0, 3);

	// 파일이 정상적으로 열렸으면
	if (file.is_open()) {
		for (int i = 0; i < MAX_NUM; i++) {
			// 0부터 3 사이의 난수 생성
			random_num = dis(gen);
			// 생성된 난수를 index로 사용하여 해당하는 뉴클레오타이드를 파일에 쓰기
			file << nucleic_sequnce[random_num];
			//if ((i + 1) % 50 == 0) file << "\n";  // 50개마다 줄바꿈 추가
		}
		cout << MAX_NUM << "개의 염기 서열이 input.txt에 저장되었습니다.";
		file.close(); //작성 후 파일 닫기
	}
	else {
		// 파일을 열지 못한 경우 에러 메시지 출력
		cout << "error";
	}

	return 0;
}