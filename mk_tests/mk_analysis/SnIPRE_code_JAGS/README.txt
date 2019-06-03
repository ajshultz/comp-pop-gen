Thanks for you interest in the SnIPRE method!  We hope to put together
an R package to run SnIPRE, but that may not be out very soon.  In the
mean time, I've attached a folder with all the necessary files you need
to play with the model.  It includes a sample script
(SnIPRE_example_script.R) and sample data set (data_example.txt), and a
sheet with brief definitions of the output.  I've tested the code and it
runs.  So, to trouble shoot I suggest you start with the included
simulated example data just to make sure it works.  Then apply to your own data.

All you should need to edit is the "SnIPRE_example_script.R".

The empirical Bayes version should be able to run on any machine, so
long has you have the "lme4" package installed in R.

I've recently changed the Bayesian implementation to utilize JAGS
instead of WinBUGS/OpenBUGS.  That way no matter what operating system
you use you do not need to have a Windows emulator.

To install JAGS, download the files from here:
http://sourceforge.net/projects/mcmc-jags/files/

Install the corresponding package for R:  R2jags  (this can be done
through the usual R interface).

Everything else you should need is in the scripts provided in the
attached folder.

The code automatically outputs estimates of the parameters in the PRF
framework (assuming neutral demography).  If you think those assumptions
are not valid then you should ignore these estimates for tau, theta,
gamma and f.   Even in this case, however, the ".class", ".est",
".lbound", ".ubound" and ".pneg" will still provide valid estimates.

Let me know if you have questions!

eilertson@psu.edu

-Kirsten Eilertson