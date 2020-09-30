#include <iostream>
#include <stdio.h>
#include <cstdlib>

int main()
{

int original_first = 419;//1;//11496;//419;//749
int original_last = 627;//209;//11550;//627;//814
int new_first = 56436;//56228;//11551;//5952;//6161

int new_last = original_last - original_first + original_last;

int ii;

char command[128];

for(ii = original_first; ii <= original_last; ii++)
{
 fprintf(stderr, "copying %i to %i\n", ii, new_first - original_first + ii);
 sprintf(command, "cp restart/restart_new_%6.6i restart/restart_new_%6.6i", ii, new_first - original_first + ii);
 system(command);
}

return 0;
}

