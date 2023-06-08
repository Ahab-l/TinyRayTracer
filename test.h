#include <stdio.h>

class Trace {
public:
	Trace()
	{
		noisy = 0;
	}
	void print(char* s); //类内没有显示声明
	void on() { noisy = 1; }
	void off() { noisy = 0; }
private:
	int noisy;
};
//类外显示定义
inline void Trace::print(char* s)
{
	if (noisy)
	{
		printf("%s", s);
	}
}
