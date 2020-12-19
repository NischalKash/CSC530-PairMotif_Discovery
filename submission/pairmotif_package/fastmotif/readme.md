# STEPS TO IMPORT THE PACKAGE

1. Move into the directory where the package is. In this case move into src folder under fastmotif directory.
2. Open the Julia Shell in the Terminal
3. Type ']' and press enter
4. You will enter the pkg shell now
5. Enter the following command : "activate ."
6. Press backspace to return to the Julia Shell
7. Enter the command : "import fastmotif"

Now you should be able to make function calls to the functions in the package


# Functions in the Fastmotif Package and its detailed description:


#### function d_hamm(x,x1)
<ul>
  <li>Calculates the hamming distance of 2 strings x, x1</li>
  <li>Input - string x, string x1</li>
  <li>Output - integer representing the hamming distance between x and x1</li>
  </ul>


#### function get_C(x, dna, d)
<ul>
<li>Calculates the set C between motif x and the rest of the dna segments with a hamming distance of d. The set C comprises of all k-mers in the rest of the segments for which the d_hamm(x, k-mer) <= 2d</li>
<li>Input - string x, list of string dna and integer d</li>
<li>Output - Dictionary Cxsi containing the set C with the dna segment as keys</li>
</ul>

#### function filter1(x, xp, z, d)
<ul>
<li>Evaluates if the strings x, xp and z conform with the First filtering rules of the paper.</li>
<li>Input - string x, string xp, string z and integer d</li>
<li>Output - true if the values satisfy the rule, false if not</li>
</ul>

#### function p_val(x,x1,z,d)
<ul>
<li>Calculates the P00 and P10 values between the strings x, xp and z. The definitions of these values can be found in the paper</li>
<li>Input - string x, string xp, string z and integer d</li>
<li>Output - List [P00,P10]</li>
</ul>

#### function filter2(x,x1,z,d)
<ul>
<li>Evaluates if the strings x, xp and z conform with the First filtering rules of the paper.</li>
<li>Input - string x, string xp, string z and integer d</li>
<li>Output - true if the values satisfy the rule, false if not</li>
</ul>

#### function get_all_md_kmers(x, xp, d, ch, k, Md)
<ul>
<li>Populates the set Md which are all k-mers that are within hamming distance bounds of x and xp and have atleast 1 y in each segment of the DNA</li>
<li>Input - string x, string xp, integer d, list ch, integer k and list Md</li>
<li>Output - Void but populates Md</li>
</ul>

#### function verify(y, Cp, r, d)
<ul>
<li>Calculates the k-mers that conform with the final rule of having atleast 1 yp in C with hamming distance <=d in each segment</li>
<li>Input - string y, dictionary Cp, integer d, integer r and integer d</li>
<li>Output - List [boolean, Dictionary]. The boolean contains whether y conforms to the verification logic and final_dict contains the position of yi in every segment</li>
</ul>

#### function pairmotif(dna, k, d_thres)
<ul>
<li>Baseline PairMotif algorithm implementation</li>
<li>Input - list of strings dna, integer k and integer d_thres</li>
<li>Output - Tuple (M, motif_positions). M is the list of consensus motifs found in the dna strands and motif_positions is a dictionary containing the actual k-mers in each segment along with the starting position of the k-mers in that segment</li>
</ul>

#### function init_screen(DNA, k, n, use_prob)
<ul>
<li>Pre-filtering on the segment S1 that runs before the main loop of PairMotif and generates the n number of initial candidate motifs. Contains implementation of 2 strategies "Entropy" and "Probability" (Both perform the same)</li>
<li>Input - list of strings dna, integer k, integer n and boolean use_prob</li>
<li>Output - List of strings which represent the initial candidate motifs from S1</li>
</ul>

#### function pairmotif_faster(dna, k, expected_n, use_prob, d_thres)
<ul>
<li>Modifier PairMotif implementation which makes use of the init_screen function above</li>
<li>Input - list of strings dna, integer k and integer d_thres</li>
<li>Output - Tuple (M, motif_positions). M is the list of consensus motifs found in the dna strands and motif_positions is a dictionary containing the actual k-mers in each segment along with the starting position of the k-mers in that segment</li>
</ul>

#### function apriori(motifs)
<ul>
<li>Post Processing function to be added to pairmotif that further prunes the motifs</li>
<li>Input - list of strings which are the consensus motifs</li>
<li>Output - Pruned list of strings/motifs which represent a shorted set of motifs</li>
</ul>

#### function get_score(motifs)
<ul>
<li>Calculate the consensus score for the motifs</li>
<li>Input - list of strings which are the motifs</li>
<li>Output - An integer value representing the consensus score of the motifs</li>
</ul>
