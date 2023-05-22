#include <stdio.h>
#include <stdlib.h>

void Abc_Start();
void Abc_Stop();

typedef struct Abc_Frame_t_ Abc_Frame_t;

Abc_Frame_t *Abc_FrameGetGlobalFrame();
int Cmd_CommandExecute(Abc_Frame_t *pAbc, const char *sCommand);

int main(int argc, char *argv[]) {
    Abc_Frame_t *pAbc;
    char *pFileName;
    char Command[1000];

    if (argc != 2) {
        printf("Wrong number of command-line arguments.");
        return -1;
    }

    pFileName = argv[1];

    Abc_Start();
    pAbc = Abc_FrameGetGlobalFrame();

    sprintf(Command, "bmatch -v -r benchmark/case%d/cir1.v benchmark/case%d/cir2.v", atoi(pFileName), atoi(pFileName));
    if (Cmd_CommandExecute(pAbc, Command)) {
        printf("Cannot execute command \"%s\".\n", Command);
        return 1;
    }

    Abc_Stop();
    return 0;
}
