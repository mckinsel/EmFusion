#include "main.h"
#include "string.h"

int main(int argc, char * argv[]){

	if (strcmp(argv[1],"EM") == 0) {
		return EM_main(argc-1, argv+1);
	}



}
