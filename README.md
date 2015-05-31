# HePhysics
Data analysis for high-energy physics (HEP). 

<h1>Introduction</h1>

HePhysics package is a Java library used for various data analysis tasks in  high-energy particle physics (HEP), which complements FreeHep and extends it for future pp collider experiments. It is a third-party module of 
<a href="http://jwork.org/dmelt/">DataMelt</a>.
In particular, it includes various Lorentz-type particles and transformations for high-energy physics applications, Java implementation of the longitudinally invariant kT,  anti-kT and Cambridge/Aachen clustering algorithms using the recombination
method.<br>
Examples are:
<ul>
<li>
KTjet - uses pseudorapidity and float values. A fast light-weight implementation
</li>

<li>
SCJet - uses rapidity for merging of particles and double values. 
</li>

<li>
SCJet - uses rapidity for merging of particles and double values.
</li>

<li>
JetN2 - Java implementation of the jet algorithm by I.Pogrebnyak using the N^2 FastJet approach.
</li>

</ul>

Stand-alone kT algorithms using the above codding, together with Java and C++ implementations, can be found in the 
<a href="https://github.com/chekanov/hepjet/">HepJet project</a>. 

<h1>Current usage</h1>

HePhysics  is a contributed package to the <a href="http://jwork.org/dmelt/">DataMelt</a> data-analysis environment and can be used in Java, Jython and Groovy.  It also used for vonline alidation analysis scripts by the <a href="http://atlaswww.hep.anl.gov/hepsim/">HepSim</a> Monte Carlo database for current and planned particle experiments.

If you have noticed bugs and  made changes, please  inform  dmelt@jwork.org .
In this case, your changes will be included to the next SCaVis release.

<h1>How to compile and test</h1>

You need to install Java 7 JDK (or above) and Athe ANT tool  to build the package. Then compile and run a simple test as:
To compile, copy the "lib" directory from the <a href="http://jwork.org/dmelt/">DataMelt</a>.
Then:

<pre>
ant
ant run
</pre>

This will run a simple test: the anti-kT jet algorithm using the input file jets/single-event.dat.
The main code in hephysics/jet/SCJet.java
The implementation is closely follows the description in  the paper "The anti-k_t jet clustering algorithm
" M.Cacciari, G.Salam,  G.Soyez, http://arxiv.org/abs/0802.1189

The input data jets/single-event.dat for this test is taken from the official  FastJet http://fastjet.fr/ web page. The C++ code is also included (in the directory "jets").

<h1>Java API</h1>
<ol>
<li><a href="http://jwork.org/dmelt/api/doc.php/hephysics/jet/package-summary">kT-type jet algorithms</a>
</li>
<li><a href="http://jwork.org/dmelt/api/doc.php/hephysics/particle/package-summary">Lorentz particles</a>
</li>
<li><a href="http://jwork.org/dmelt/api/doc.php/hephysics/hepsim/PromcUtil">package to fill  particle array from ProMC format</a>
</li>
</ol>

