#pragma once

#ifdef SMITH_INFO
 #define I_
#else
 #define I_ if (0)
#endif

#ifdef SMITH_DEBUG
 #define D_
#else
 #define D_ if (0)
#endif