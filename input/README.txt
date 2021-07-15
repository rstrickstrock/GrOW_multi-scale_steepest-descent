Installation and execution of GrOW:

1. copy the "Code" directory to the desired location
2. navigate to Code/GrOW/
3. exectute the main script the following way:
	python main.py <config>
	e.g.
	$ python main.py ../opt_1/octane_hybrid.cfg



Tips:

- make sure all the relative filepaths in the <config>-file are correct
- make sure that all necessary pyhton packages are installed
- adapt files in Code/GrOW/parallel_jobs/*
		 Code/GrOW/simulation/*
  	to fit to your simulation tools and environments used (e.g. cluster queueuing software)


