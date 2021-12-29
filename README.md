ntreeshape C++ code to calculate the number of trees of a given size
From the number of trees, the number of genetic programming programs
of that size can be calculated.

Distribution shown here 
http://www.cs.ucl.ac.uk/staff/W.Langdon/fitnessspaces/node3.html
http://www.cs.ucl.ac.uk/staff/W.Langdon/fitnessspaces/node4.html

Code used in 
A Field Guide to Genetic Programming, ISBN 978-1-4092-0073-4
http://www.gp-field-guide.org.uk/


<p>
<hr>
<P>
Contents
<UL>
<li>ntrees.cc, Apr  9  2001, r1.2
<br>
Calculate number of trees of a given size

<li>ntreeshape.cc, Aug 20  2002, r1.22
<br>Calculate number of trees by depth for a given size

<li><a href="README.md">README.md</a>
<br>This file

</ul>

<!--hmm that was shit
 \pi
$\pi$-->
&pi;
<span>&#8730;</span>

Note the mean and variance for various types of tree 
have been given asymptotically for large trees
by Flajolet and co-authors.
E.g. for random binary trees:

mean depth = $<span>&#8730;2&pi;(n-1)</span>

variance = 4&pi;(&pi;/3-1)(n-1)/2

Note the formulae are often given in terms of
number of internal nodes, rather than tree size (n, including leaf nodes).
And for trees under 10000 nodes, the actual mean
and variance can be somewhat different from the 
asymptotic formula.



As an example of software bit-rot, 
consider the differences between the first versions
of these C++ files in GitHub (dating from 2001 and 2002)
and the current versions.
http://blog.ieeesoftware.org/2020/03/bit-rot-computer-software-degrades-over.html
