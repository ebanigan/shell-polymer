/* This subroutine runs if MEASURE_T_ATTACH is true.  It prints out the values in the vectors t_att_bins and t_strong_att_bins.
 * This tells us the frequency of different attachment and "strong" (U > kT) attachment times in units of time steps.
 */
#if TETHER
#if MEASURE_T_ATTACH
void print_times()
{
long ii, jj;

FILE *attachfile;
FILE *strongattachfile;
char attachname[72];
char strongattachname[72];
//FILE *num_s_attachments;
//char numname[72];
FILE *detachfile;
char detachname[72];
FILE *weakdetachfile;
char weakname[72];

sprintf(attachname, "output/time_attached%6.6li", TRIALNUMBER);
sprintf(strongattachname, "output/time_strong_attached%6.6li", TRIALNUMBER);
attachfile = fopen(attachname, "w");
strongattachfile = fopen(strongattachname, "w");

for(int qrst = 0; qrst < t_att_bins.size(); qrst++)
  fprintf(attachfile, "%li %li\n", t_att_bins[qrst], t_att_count[qrst]);
for(int lmnop = 0; lmnop < t_strong_att_bins.size(); lmnop++)
  fprintf(strongattachfile, "%li %li\n", t_strong_att_bins[lmnop], t_strong_att_count[lmnop]);

for(int abc1 = 0; abc1 < NUMBER_OF_MONOMERS; abc1++)
if((mono_list[abc1]).get_t_attached() > 0)
{
fprintf(stderr, "%li %li\n", (mono_list[abc1]).get_t_attached(), (mono_list[abc1]).get_t_strong_attached());
}


fflush(attachfile);
fclose(attachfile);
fflush(strongattachfile);
fclose(strongattachfile);

sprintf(detachname, "output/time_detached%6.6li", TRIALNUMBER);
sprintf(weakname, "output/time_weak_detached%6.6li", TRIALNUMBER);

weakdetachfile = fopen(weakname, "w");
detachfile = fopen(detachname, "w");

for(ii = 0; ii < t_detach_bins.size(); ii++)
  fprintf(detachfile, "%li %li\n", t_detach_bins[ii], t_detach_count[ii]);
for(ii = 0; ii < t_weak_detach_bins.size(); ii++)
  fprintf(weakdetachfile, "%li %li\n", t_weak_detach_bins[ii], t_weak_detach_count[ii]);

fprintf(stderr, "detach times:\n");
for(ii = 0; ii < NUMBER_OF_MONOMERS; ii++)
if((mono_list[ii]).get_t_weak_detached() > 0)
{
fprintf(stderr, "%li %li\n", (mono_list[ii]).get_t_detached(), (mono_list[ii]).get_t_weak_detached());
}


fflush(detachfile);
fclose(detachfile);
fflush(weakdetachfile);
fclose(weakdetachfile);

/*sprintf(numname, "output/num_s_attachments%6.6li", TRIALNUMBER);
num_s_attachments = fopen(numname, "w");
for(ii = 0; ii < NUMBER_OF_MONOMERS; ii++)
{
  for(jj = 0; jj < (mono_list[ii]).num_strong_att_periods_size(); jj++)
     fprintf(num_s_atachments, "%li ", (mono_list[ii]).get_num_strong_att_periods(jj));}
  fprintf(num_s_atachments, "\n");
}
fflush(num_s_attachments);
fclose(num_s_attachments);*/
}
#endif
#endif
