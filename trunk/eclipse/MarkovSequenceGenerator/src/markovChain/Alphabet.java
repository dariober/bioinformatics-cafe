package markovChain;

import java.util.HashSet;
import java.util.Set;

/**
 * Class for the various alphabets that might be supported. 
 * E.g. nucleotides A C G T or A C T G a c t g, or
 * aminoacids, etc.
 * @author berald01
 *
 */
public class Alphabet {

	private Set<Character> nucleotides= new HashSet<Character>();
	
	public Set<Character> getNucleotides() {

		this.nucleotides.add('A');
		this.nucleotides.add('C');
		this.nucleotides.add('G');
		this.nucleotides.add('T');
		
		return nucleotides;
	}
	
	/**
	 * Produce a set of unique elements in input sequence that will be
	 * used as alphabet by the Markov model.
	 * @param sequence
	 * @return
	 */
	public Set<Character> getAlphabetFromString(String sequence){
		Set<Character> alphabet= new HashSet<Character>();
		for (char x : sequence.toCharArray()){
			alphabet.add(x);
		}
		return(alphabet);
	}
}


