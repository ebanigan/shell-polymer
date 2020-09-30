
vector<double> quicksort(vector<double> xlist)
{
int ii;

double x0;
vector<double> lesslist;
vector<double> greaterlist;
x0 = xlist[0];
vector<double> outlist;

for(ii = 1; ii < xlist.size(); ii++)
{
  if(xlist[ii] < x0)//+1e-7*(ranf0()-0.5))//for this code, shouldn't need to worry about identical entries.
    lesslist.push_back(xlist[ii]);
  else
    greaterlist.push_back(xlist[ii]);
}

if(lesslist.size() > 1)
 lesslist = quicksort(lesslist);
if(greaterlist.size() > 1)
 greaterlist = quicksort(greaterlist);

for(ii = 0; ii < lesslist.size(); ii++)
 outlist.push_back(lesslist[ii]);

outlist.push_back(x0);

for(ii = 0; ii < greaterlist.size(); ii++)
 outlist.push_back(greaterlist[ii]);

return(outlist);
}




