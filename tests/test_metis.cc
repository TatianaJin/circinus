// Copyright 2021 HDL
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
#include "gflags/gflags.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

#include "metis.h"
#include <stdio.h>
/*
 0--1--2
 |  |  |
 3--4--5
 |  |  |
 6--7--8
*/
int main()
//TEST(TestMetis, SimpleGraph)
{
    idx_t xadj[]={0,2,5,7,10,14,17,19,22,24};
    idx_t adjncy[]={1,3,0,2,4,1,5,0,4,6,1,3,5,7,2,4,8,3,7,4,6,8,5,7};
    idx_t nvtxs=9;
    idx_t ncon=1;
    idx_t nparts=4;
    idx_t objval;
    idx_t part[9];
    printf("%d\n",METIS_PartGraphKway(&nvtxs,&ncon,xadj,adjncy,NULL,NULL,NULL,&nparts,NULL,NULL,NULL,&objval,part));
    printf("cut-edge count=%d\n",objval);
    for(int i=0;i<9;++i)printf("%d ",part[i]);
    printf("\n");
}
