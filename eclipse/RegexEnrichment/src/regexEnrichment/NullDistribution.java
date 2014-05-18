package regexEnrichment;

import java.util.HashMap;
import java.util.regex.Pattern;

import com.google.common.collect.TreeMultiset;

public class NullDistribution {
	
	private markovChain.Model markovModel;
	private int orderMarkovChain;
	
	/*                         C O N S T R U C T O R S               
	 * ---------------------------------------------------------------------- */
	
	NullDistribution(){
	}
	
	NullDistribution(String sequence, int orderMarkovChain, double pctPseudocounts){
		this.orderMarkovChain= orderMarkovChain;
		if (orderMarkovChain > 0){
			markovChain.Model markovModel= Tools.generateMarkovModel(
					orderMarkovChain, sequence, pctPseudocounts);
			this.markovModel= markovModel;
		}
		else {
			this.markovModel= new markovChain.Model();
		}
	}
	
	/*                              M E T H O D S                   
	 * ---------------------------------------------------------------------- */
	/**
	 * Return the nucleotide frequency of each nucleotide in input sequence. 
	 * @param sequence
	 * @param precision
	 * @return An Object array where element 0 is an array of nucleotides and 
	 * element 1 is array of frequencies. Arrays 0 and 1 align side by side.
	 */
	private Object[] pairNucleotidesAndProbs(String sequence, int precision){
		HashMap<String, Float> ntProbs= Tools.nucleotideFrequencies(sequence, precision);		
		// Convert key/values (nucleotide/frequency) to a pair of arrays for speed
		// later when getting cumulative prob. func.
		String[] nucleotides= new String[ntProbs.size()];
		float[] nucleotProb= new float[ntProbs.size()];
		int x= 0;
		for (String k : ntProbs.keySet()){
			nucleotides[x]= k;
			nucleotProb[x]= ntProbs.get(k);
			x++;
		}
		return new Object[] {nucleotides, nucleotProb};
	}
	
	/** 
	 * Shuffle the input string and search for the given pattern to generate a null distribution
	 * of hits produced by chance.
	 * @param sequence Input sequence to shuffle
	 * @param patternString The regex pattern to search in the shuffled sequences 
	 * @param nrnd Number of randomizations to execute
	 * @param precision see above 
	 * @param Counter of number of hits. Possibly from  previous loop which needs to be incremented  
	 * @param revcomp Should the input sequence be reverse complemented to search for pattern?  
	 * @return Dictionary where the key is the number of hits within the same sequence and the value
	 * 		is the number of sequences found. E.g. [0 x 19069, 1 x 927, 2 x 4] means: 
	 *      19069 sequences had 0 regex hits; 927 seqs had exactly 1 hit; 4 seqs had 2 hits.
	 */
	public TreeMultiset<Integer> generateNullDistribution(
			String sequence, Pattern pattern, int nrnd, int precision, int obsHits,
			TreeMultiset<Integer> nMatches, boolean revcomp){
		
		int i= 1 + nMatches.size();
		float pvalObsHits= -1;
		while(nrnd >= i){
			// Generate ranodm sequence according to model:
			String xsequence;
			if (orderMarkovChain > 0){
				xsequence= new markovChain.MarkovSequence().generateMarkovSequence(
						markovModel, sequence.length());
			} 
			else {
				Object[] pairArray= pairNucleotidesAndProbs(sequence, precision);
				String[] nucleotides= (String[])pairArray[0];
				float[] nucleotProb= (float[])pairArray[1];
				xsequence= Tools.generateRandomSequence(nucleotides, nucleotProb, sequence.length());
			}
			int nm= Tools.regexCounter(pattern, xsequence, revcomp);
			nMatches.add(nm);
			// Every n iterations check if the pvalue is greater then threshold. 
			// If so, there is no point in going any further, the observed hits are
			// non significant.
			if( i % 500 == 0 && i >= 1000){
				pvalObsHits= Tools.mapObservedHitsToNull(nMatches, obsHits);
				if (pvalObsHits > 100 / (float) i) {
					break;
				}
			}
			i++;
		}
		return(nMatches);
	}

	public String toString(){
		StringBuilder sb= new StringBuilder();
		sb.append("Order as set: " + this.orderMarkovChain + "\n");
		try{
			sb.append("Markov model:\n" + this.markovModel.toString());
		}
		catch(NullPointerException e){
			sb.append("Markov model: " + "N/A");
		}
		return(sb.toString());
	}
	
	/*                       S E T T E R S   &   G E T T E R S                   
	 * ---------------------------------------------------------------------- */
	
	public markovChain.Model getMarkovModel() {
		return markovModel;
	}

	public void setMarkovModel(markovChain.Model markovModel) {
		this.markovModel = markovModel;
	}

	public int getOrderMarkovChain() {
		return orderMarkovChain;
	}

	public void setOrderMarkovChain(int orderMarkovChain) {
		this.orderMarkovChain = orderMarkovChain;
	}
	
}
