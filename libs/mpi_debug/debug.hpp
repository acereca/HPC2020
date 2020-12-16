#include <iostream>
#include <unistd.h>

namespace MPI_DBG {
    /*
    *  enables debugging with gdb
    *  (after https://www.open-mpi.org/faq/?category=debugging#serial-debuggers)
    *
    *  - requires MPI_Init beforehand
    *  - to connect run:
    *    $ gdb attach <pid>
    *    (gdb) break sleep
    *    (gdb) frame 3
    *    (gdb) set var ifl = 1
    *    (gdb) break ...
    */
    inline void enable_debugging(int rank) {
        #ifdef DBG
            char host[256];
            gethostname(host, sizeof(host));
            std::cout << "Rank " << rank << ": Host " << host << ", PID " << getpid() << ", ready to attach!" << std::endl;
            volatile int ifl = 0;
            while (0 == ifl) sleep(5);
        #endif
    }
}
