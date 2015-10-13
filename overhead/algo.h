/*
Copyright (c) Durham University, 2013

This file is part of PARMES.

PARMES is free software; you can redistribute it
and/or modify it under the terms of the GNU General Public License
version 3 as published by the Free Software Foundation.

PARMES is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PARMES. If not, see <http://www.gnu.org/licenses/>.
*/

#define DRAND() ((REAL) rand () / (REAL) RAND_MAX)

#define MUL(a, eps, b)\
{\
  (b) [0] = (a) [0] * (eps);\
  (b) [1] = (a) [1] * (eps);\
  (b) [2] = (a) [2] * (eps);\
}

#define MUL4(a, eps, b)\
{\
(b) [0] = (a) [0] * (eps);\
(b) [1] = (a) [1] * (eps);\
(b) [2] = (a) [2] * (eps);\
(b) [3] = (a) [3] * (eps);\
}

#define MUL44(a, eps, b)\
{\
(b) [0] = (a) [0] * (eps);\
(b) [1] = (a) [1] * (eps);\
(b) [2] = (a) [2] * (eps);\
(b) [3] = (a) [3] * (eps);\
(b) [4] = (a) [4] * (eps);\
(b) [5] = (a) [5] * (eps);\
(b) [6] = (a) [6] * (eps);\
(b) [7] = (a) [7] * (eps);\
(b) [8] = (a) [8] * (eps);\
(b) [9] = (a) [9] * (eps);\
(b) [10] = (a) [10] * (eps);\
(b) [11] = (a) [11] * (eps);\
(b) [12] = (a) [12] * (eps);\
(b) [13] = (a) [13] * (eps);\
(b) [14] = (a) [14] * (eps);\
(b) [15] = (a) [15] * (eps);\
}

#define DOT(a, b)\
  ((a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2])

#define ADD44(a, b, c)\
{\
    (c) [0] = (a) [0] + (b) [0];\
    (c) [1] = (a) [1] + (b) [1];\
    (c) [2] = (a) [2] + (b) [2];\
    (c) [3] = (a) [3] + (b) [3];\
    (c) [4] = (a) [4] + (b) [4];\
    (c) [5] = (a) [5] + (b) [5];\
    (c) [6] = (a) [6] + (b) [6];\
    (c) [7] = (a) [7] + (b) [7];\
    (c) [8] = (a) [8] + (b) [8];\
    (c) [9] = (a) [9] + (b) [9];\
    (c) [10] = (a) [10] + (b) [10];\
    (c) [11] = (a) [11] + (b) [11];\
    (c) [12] = (a) [12] + (b) [12];\
    (c) [13] = (a) [13] + (b) [13];\
    (c) [14] = (a) [14] + (b) [14];\
    (c) [15] = (a) [15] + (b) [15];\
}

#define ADD4(a, b, c)\
{\
  (c) [0] = (a) [0] + (b) [0];\
  (c) [1] = (a) [1] + (b) [1];\
  (c) [2] = (a) [2] + (b) [2];\
  (c) [3] = (a) [3] + (b) [3];\
}

#define ADD(a, b, c)\
{\
  (c) [0] = (a) [0] + (b) [0];\
  (c) [1] = (a) [1] + (b) [1];\
  (c) [2] = (a) [2] + (b) [2];\
}

#define SUB(a, b, c)\
{\
  (c) [0] = (a) [0] - (b) [0];\
  (c) [1] = (a) [1] - (b) [1];\
  (c) [2] = (a) [2] - (b) [2];\
}


#define SUB4(a, b, c)\
{\
  (c) [0] = (a) [0] - (b) [0];\
  (c) [1] = (a) [1] - (b) [1];\
  (c) [2] = (a) [2] - (b) [2];\
  (c) [3] = (a) [3] - (b) [3];\
}

