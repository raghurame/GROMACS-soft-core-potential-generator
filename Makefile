all:
	gcc -o generateSoftPotentialTable generateSoftPotentialTable.c -lm -frounding-math -fsignaling-nans;
	rm -rf table;
	mkdir table;
	mkdir table/soft0.01 table/soft0.05 table/soft0.08 table/soft0.1 table/soft0.2 table/soft0.3 table/soft0.4 table/soft0.5 table/soft0.6 table/soft0.7 table/soft0.8 table/soft0.9 table/soft1.0;
	./generateSoftPotentialTable table/soft0.01/table.xvg 2 0.01 0.5 7.80E-03 3.24314E-05;
	./generateSoftPotentialTable table/soft0.05/table.xvg 2 0.05 0.5 7.80E-03 3.24314E-05;
	./generateSoftPotentialTable table/soft0.08/table.xvg 2 0.08 0.5 7.80E-03 3.24314E-05;
	./generateSoftPotentialTable table/soft0.1/table.xvg 2 0.1 0.5 7.80E-03 3.24314E-05;
	./generateSoftPotentialTable table/soft0.2/table.xvg 2 0.2 0.5 7.80E-03 3.24314E-05;
	./generateSoftPotentialTable table/soft0.3/table.xvg 2 0.3 0.5 7.80E-03 3.24314E-05;
	./generateSoftPotentialTable table/soft0.4/table.xvg 2 0.4 0.5 7.80E-03 3.24314E-05;
	./generateSoftPotentialTable table/soft0.5/table.xvg 2 0.5 0.5 7.80E-03 3.24314E-05;
	./generateSoftPotentialTable table/soft0.6/table.xvg 2 0.6 0.5 7.80E-03 3.24314E-05;
	./generateSoftPotentialTable table/soft0.7/table.xvg 2 0.7 0.5 7.80E-03 3.24314E-05;
	./generateSoftPotentialTable table/soft0.8/table.xvg 2 0.8 0.5 7.80E-03 3.24314E-05;
	./generateSoftPotentialTable table/soft0.9/table.xvg 2 0.9 0.5 7.80E-03 3.24314E-05;
	./generateSoftPotentialTable table/soft1.0/table.xvg 2 1.0 0.5 7.80E-03 3.24314E-05;

