#include <iostream>
#include <climits>

using namespace std;

#define MSIZE 128

int main(){
  int i;
  int line = 0;
  int minimum[4] = {INT_MAX, INT_MAX, INT_MAX, INT_MAX};
  int maximum[4] = {0,0,0,0};
  int sum[4] = {0,0,0,0};
  int in[4];
  int a,b,c,d;
  while(cin >> a >> b >> c >> d){
    line++;

    in[0] = a; in[1] = b;
    in[2] = c; in[3] = d;

    for(i=0; i<4; i++){
      if(minimum[i] > in[i])
        minimum[i] = in[i];
      if(maximum[i] < in[i])
        maximum[i] = in[i];
      sum[i] += in[i];
    }
  }

  for(i=0; i<4; i++){
    double avg = static_cast<double>(sum[i])/static_cast<double>(line);
    cout << MSIZE*MSIZE - maximum[i] << " " << MSIZE*MSIZE - minimum[i] << " " << MSIZE*MSIZE - avg << endl;
  }

	return 0;
}
