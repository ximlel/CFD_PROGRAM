#!/bin/bash

ulimit -s  102400

cd ..
make
cd ./src

read

#./hydrocode.out Sod_10_test2 Sod_10_test2/Sod_10_test_ROE 2 1_Roe t1

#./hydrocode.out Sod_test Sod_test/Sod_test_ROE	1 1_Roe		  free
#./hydrocode.out Sod_test Sod_test/Sod_test	1 1_Riemann_exact free
#./hydrocode.out Sod_test Sod_test/Sod_test	1 2_GRP		  free
#./hydrocode.out Sod_10_test Sod_10_test/Sod_10_test_ROE 2 1_Roe		  Sod
#./hydrocode.out Sod_10_test Sod_10_test/Sod_10_test	 2 1_Riemann_exact Sod
#./hydrocode.out Sod_10_test Sod_10_test/Sod_10_test	 2 2_GRP		  Sod

#./hydrocode.out odd_even odd_even/odd_even_Roe 2 1_Roe 	  odd_even
./hydrocode.out odd_even odd_even/odd_even	2 1_Riemann_exact Sod odd_even

#./hydrocode.out NEW_TEST_OBLIQUE NEW_TEST_OBLIQUE/NEW_TEST_OBLIQUE_Roe	2 1_Roe			oblique_periodic
#./hydrocode.out NEW_TEST_OBLIQUE NEW_TEST_OBLIQUE/NEW_TEST_OBLIQUE	2 1_Riemann_exact	oblique_periodic
#./hydrocode.out NEW_TEST_OBLIQUE NEW_TEST_OBLIQUE/NEW_TEST_OBLIQUE	2 2_GRP			oblique_periodic
#./hydrocode.out NEW_TEST_OBLIQUE_SMALL NEW_TEST_OBLIQUE_SMALL/NEW_TEST_OBLIQUE_SMALL	2 1_Riemann_exact			oblique_periodic

#./hydrocode.out NEW_TEST NEW_TEST/NEW_TEST_Roe	2 1_Roe			free
#./hydrocode.out NEW_TEST NEW_TEST/NEW_TEST	2 1_Riemann_exact	free
#./hydrocode.out NEW_TEST_BIG NEW_TEST_BIG/NEW_TEST_BIG_Roe	2 1_Roe			free

#./hydrocode.out Riemann_2D3_Quad Riemann_2D3_Quad/Riemann_2D3_Quad_Roe	2 1_Roe			free
#./hydrocode.out Riemann_2D3_Quad Riemann_2D3_Quad/Riemann_2D3_Quad	2 1_Riemann_exact	free
#./hydrocode.out Riemann_2D3_Quad Riemann_2D3_Quad/Riemann_2D3_Quad	2 2_GRP			free

#./hydrocode.out RMI/RMI_81	RMI/RMI_81	2 1_Riemann_exact RMI
#./hydrocode.out RMI/RMI_81	RMI/RMI_81	2 2_GRP		  RMI
#./hydrocode.out RMI/RMI_321	RMI/RMI_321	2 1_Riemann_exact RMI
#./hydrocode.out RMI/RMI_321	RMI/RMI_321	2 2_GRP		  RMI
