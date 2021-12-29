ntreeshape C++ code to calculate the number of trees of a given size
From the number of trees, the number of genetic programming programs
of that size can be calculated.

Distribution shown here 
http://www.cs.ucl.ac.uk/staff/W.Langdon/fitnessspaces/node3.html
http://www.cs.ucl.ac.uk/staff/W.Langdon/fitnessspaces/node4.html

Code used in 
"Scaling of Program Tree Fitness Spaces", W.B.Langdon, Evolutionary Computation 7(4) pp399-428 doi:10.1162/evco.1999.7.4.399
and
"A Field Guide to Genetic Programming", ISBN 978-1-4092-0073-4
http://www.gp-field-guide.org.uk/


<p>
<hr>
<P>
Contents
<UL>
<li>ntrees.cc, Apr  9  2001, r1.2
<br>
Calculate number of trees of a given size

<li>ntrees_example141.out
<br>Expected output of ./ntrees 1,141,2 1 0 1
<br>Note by having only one leaf and one binary function
the number of programs is the same as the number of tree shapes

<li>ntreeshape.cc, Aug 20  2002, r1.22 (now r1.23)
<br>Calculate number of trees by depth for a given size

<li>ntreeshape_example141.out
<br>Expected output of ./ntreeshape 1,141,2 1,100

<li><a href="README.md">README.md</a>
<br>This file

</ul>

<!--dont work  \pi $\pi$-->
<!--pi and square root ok as &pi; <span>&#8730;</span> !-->
<!--Dec 2021 MathML not supported by GitHub 
https://github.com/github/markup/issues/551
<math xmlns="http://www.w3.org/1998/Math/MathML" display="block">
<mrow><msup><mrow><mi>e</mi></mrow><mrow><msqrt><mi>x</mi></msqrt></mrow></msup><mo>-</mo><mfrac><mrow><mfrac><mrow><msup><mrow><mo>sin</mo></mrow><mrow><mo>-</mo><mn>1</mn></mrow></msup><mspace width="0.167em"></mspace><mfenced open="(" close=")" separators=""><mn>2</mn><mspace width="0.167em"></mspace><mi>x</mi></mfenced></mrow><mrow><mn>2</mn><mo>&times;</mo><msup><mrow><mn>10</mn></mrow><mrow><mn>10</mn></mrow></msup><mo>+</mo><msup><mrow><mi>x</mi></mrow><mrow><mn>3</mn></mrow></msup></mrow></mfrac></mrow><mrow><mo>-</mo><mn>12</mn></mrow></mfrac></mrow>
</math>

<msqrt> base </msqrt>
<mroot> base index </mroot>

!-->

Note the mean and variance for various types of tree 
have been given asymptotically for large trees
by Flajolet and co-authors.
E.g. for random binary trees:

mean depth = (2&pi;(n-1))<sup>0.5</sup>

variance = 4&pi;(&pi;/3-1)(n-1)/2

Notice for random binary trees, although the distribution tends
to a Gaussian as the trees get bigger,
it continues to have a broad spread,
with the standard deviation being 22% <!--0.21725!-->
of the mean.


Formulae are often given in terms of
number of internal nodes (N), rather than tree size
(n, including leaf nodes, for binary trees N=(n-1)/2).
Also for trees under 10000 nodes, the actual mean
and variance can be somewhat different from the 
asymptotic formula.

Note: in this version of ntreeshape.cc
depth starts at 1 (rather than zero).

As an example of software bit-rot, 
consider the differences between the first versions
of these C++ files in GitHub (dating from 2001 and 2002)
and the current versions.
http://blog.ieeesoftware.org/2020/03/bit-rot-computer-software-degrades-over.html
