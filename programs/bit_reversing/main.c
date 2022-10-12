#include <stdio.h>
#include <stdlib.h>

int reverse(int num) {

    unsigned int count = 9;
    unsigned int tmp = num;         //  Assign num to the tmp

    num >>= 1; // shift num because LSB already assigned to tmp

    while(num)
    {
       tmp <<= 1;  //shift the tmp because alread have the LSB of num

       tmp |= num & 1; // putting the set bits of num

       num >>= 1;

       count--;
    }

    tmp <<= count; //when num become zero shift tmp from the remaining counts
    tmp &= 0x000003FF;
    return tmp;
}

int main()
{
    int tomb [1024];

    for(int i=0; i<1024; i++){
        tomb[i] = reverse(i);
        printf("Tomb: ind = %d, ert=%d\n",i , tomb[i]);
    }
    return 0;
}
