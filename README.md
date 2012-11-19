<h3>About</h3>
This is the official repository for the Find RNA 3D (FR3D) program for searching
RNA 3D structures.

Originally the code was hosted at the
<a href="http://rna.bgsu.edu/FR3D">FR3D project homepage</a>,
but the development will now be moved to Github. The first 4 commits marked with
dates were taken from the project's website and were authored by
Professor <a href="http://www-math.bgsu.edu/z/">Craig Zirbel</a> and other
members of the <a href="http://rna.bgsu.edu">BGSU bioinformatics group</a>.
The transition to Github is managed by
<a href="https://github.com/AntonPetrov">Anton Petrov</a>.

<h3>References</h3>
The <a href="http://www.ncbi.nlm.nih.gov/pubmed/17694311">original publication</a>
describing FR3D and its methodology appeared in the Journal of Mathematical Biology
in 2008. In 2011 the online version of FR3D called
<a href="http://rna.bgsu.edu/webfr3d">WebFR3D</a> was made available.

<h3>Requirements</h3>
* Matlab (tested on R2007b and later). Some FR3D functionality is also available
when running in <a href="http://www.gnu.org/software/octave/">Octave</a>,
but it is not well-tested.

* Internet connection may be required for downloading 3D structures from the
<a href="http://pdb.org">PDB</a>.


<h3>Installation</h3>

1. Download the source code.

        git clone https://github.com/BGSU-RNA/FR3D.git

2. Launch Matlab and set path:

        cd /path/to/your/local/FR3D
        setup_path

3. To run FR3D GUI in Matlab type:

        FR3D


<h3>Contact Us</h3>
Feel free to <a href="https://github.com/BGSU-RNA/FR3D/issues">submit an issue</a>
or get in touch via a contact form on the website
of the <a href="http://rna.bgsu.edu">BGSU RNA Bioinformatics group</a>.