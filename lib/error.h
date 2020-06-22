#ifndef __Error__
#define __Error__

#include <vector>
#include <string>


double average(std::vector<double> v);
double error(std::vector<double> v);

void do_blocking_method(std::string file_in, std::string file_out);
void do_raw_blocking_method(std::string file_in, std::string file_out, int n_blocks);

#endif
