#ifndef __externC_H_INCLUDED
#define __externC_H_INCLUDED

extern "C" void getimdata();
extern "C" void fillinput(input_h *);
extern "C" void write_fits(flux,input_h,long,const char *);
extern "C" void exit_with_usage(int);
extern "C" void print_warning(int);

extern "C" void print_header(const char *,const char *);
extern "C" void printstat(long);
extern "C" void printdone();
extern "C" void printmem(size_t,int);
extern "C" void printmemreq(size_t,int);
extern "C" void printmemreqcuda(size_t,int,int);
extern "C" void printcudamem(size_t,int,int);
extern "C" void printinitcuda(int);
extern "C" void printinitheader(int);

extern "C" size_t getMemorySize();

extern "C" int read_in_param(char *);
extern "C" int nextwrite(int);

#endif
                                                          