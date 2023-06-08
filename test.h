#include <stdio.h>

class Trace {
public:
	Trace()
	{
		noisy = 0;
	}
	void print(char* s); //����û����ʾ����
	void on() { noisy = 1; }
	void off() { noisy = 0; }
private:
	int noisy;
};
//������ʾ����
inline void Trace::print(char* s)
{
	if (noisy)
	{
		printf("%s", s);
	}
}
