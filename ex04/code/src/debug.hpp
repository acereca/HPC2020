#include <iostream>
#include <unistd.h>

namespace MPI_DBG {
    inline void enable_debugging(int rank) {
        #ifdef DBG
            std::cout << "Rank " << rank << ": PID " << getpid() << ", ready to attach!" << std::endl;
            volatile int ifl = 0;
            while (0 == ifl) sleep(5);
        #endif
    }
}
