#include <stdio.h>
#include <stdlib.h>

void Abc_Start();
void Abc_Stop();

typedef struct Abc_Frame_t_ Abc_Frame_t;

Abc_Frame_t *Abc_FrameGetGlobalFrame();
int Cmd_CommandExecute(Abc_Frame_t *pAbc, const char *sCommand);

int main(int argc, char *argv[]) {
    Abc_Frame_t *pAbc;
    char Command[1000];

    if (argc != 3) {
        printf("Usage: %s <input> <match>\n", argv[0]);
        return -1;
    }

    Abc_Start();
    pAbc = Abc_FrameGetGlobalFrame();

    sprintf(Command, "bmatchgroup -v -r %s %s", argv[1], argv[2]);
    if (Cmd_CommandExecute(pAbc, Command)) {
        printf("Cannot execute command \"%s\".\n", Command);
        return 1;
    }

    Abc_Stop();
    return 0;
}
