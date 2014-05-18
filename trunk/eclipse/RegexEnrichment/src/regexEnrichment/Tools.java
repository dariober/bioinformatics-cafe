package regexEnrichment;

import java.io.*;
import java.lang.StringBuilder;

import com.google.common.collect.*;

import java.util.*;
import java.util.regex.*;

import org.biojava3.core.sequence.*;
import org.biojava3.core.sequence.compound.*;

public class Tools {

	private static AmbiguityDNACompoundSet ambiguityCompoundSet= new AmbiguityDNACompoundSet();
	
	/**
	 * Nucleotide composition of the input string.
	 * @param sequence String to get count of individual chars.
	 * @return Count of chars.
	 */
	public static HashMultiset<String> nucleotideComposition(String sequence){

		List<String> xlist= new ArrayList<String>();
		
		for( char c : sequence.toCharArray()){
			xlist.add(Character.toString(c));
		}
		
		HashMultiset<String> ntCount= HashMultiset.create();

		ntCount.addAll(xlist);	
		return(ntCount);
	}

	public static float roundProb(float d, int precision){
		if (precision > 0){
			float rounded= Math.round(d * (float) precision) / (float) precision;
			return(rounded);
		}
		else{
			return(d);
		}
	}
	
	/**
	 * Transform raw nucleotide counts to frequencies rounded to given precision.
	 * 
	 * If any of the nucleotides has been rounded to 0, use the exact
	 * frequencies. E.g. 
	 * sequence= 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAG'
	 * ntCount= {A= 49, G= 1}
	 * freq= {A= 0.98, G= 0.02}
	 * roundedProbWithPrecision20= {A= 1, G= 0}
	 * 
	 * The rounded probabilities will generate random sequence composed
	 * by 'A' homopolymers!
	 *
	 * NOTE: Sum of rounded % will not always sum to 1. See also:
	 * http://stackoverflow.com/questions/5227215/how-to-deal-with-the-sum-of-rounded-percentage-not-being-100
	 *
	 * 
	 * @param ntCount Count of each nucleotide.
	 * @param precision @see nucleotideFrequencies
	 * @return Rounded frequency of each nucleotide.
	 */
	public static HashMap<String, Float> countsToProbs(HashMultiset<String> ntCount, int precision){
		HashMap<String, Float> probs= new HashMap<String, Float>();
		float totCount= (float) ntCount.size();

		boolean hasZero= false;
		for( String sequence : ntCount.elementSet() ){
			float prob= roundProb(ntCount.count(sequence) / totCount, precision);
			if (prob == 0){
				// If a rounded prob gets down to zero, exit and use "exact" frequencies.
				hasZero= true;
				break;
			}
			probs.put(sequence, prob);
		}
		if (hasZero){
			for( String sequence : ntCount.elementSet() ){
				float prob= ntCount.count(sequence) / totCount;
				probs.put(sequence, prob);
			}
		}
		return(probs);
	}
	
	/**
	 * Return the nucleotide frequencies in "sequence" rounded to "precision" 
	 * @param sequence
	 * @param precision Round the nucleotide frequencies found in "sequence" to nearest 1/precision.
	 * 		E.g. precision= 20 will round 0.11 to 0.10 and 0.16 to 0.15. 
	 * @return
	 */
	public static HashMap<String, Float> nucleotideFrequencies(String sequence, int precision){
		HashMultiset<String> ntCount= Tools.nucleotideComposition(sequence);
		HashMap<String, Float> ntFreq= Tools.countsToProbs(ntCount, precision);
		return(ntFreq);
	}
	
	/**
	 * Generate a random sequence of length n using the nucleotide probabilities ntProbs.
	 * For speed consider using the other version.
	 * @param ntProbs
	 * @param n Lenght of random sequence
	 * @return
	 */
	private static Random rnd= new Random();
	public static String generateRandomSequence(HashMap<String, Float> ntProbs, int n){
		StringBuilder shuffledSeq= new StringBuilder();
		int i= 1;
		// Random rnd= new Random();
		while(i <= n){
			float r = rnd.nextFloat();
			float cdf = (float) 0.0;
			for (String k : ntProbs.keySet()) {
				cdf += ntProbs.get(k);
			    if (r <= cdf) {
			        shuffledSeq.append(k);
			        break;
			    }
			}
			i++;
		}
		return(shuffledSeq.toString());
	}

	/**
	 * Generate random sequence by sampling from nucleotides array with
	 * probability nucleotFreqs. 
	 * MEMO: arrays nucleotides and nucleotFreqs must align side by side
	 * to get the correct assignments.
	 * A HashMap might be more appropriate but it is slower. 
	 * @param nucleotides Alphabet of nucleotides to sample
	 * @param nucleotProb Probability of sampling each nucleotide
	 * @param n Lenght of random string to return
	 * @return Random sequence
	 */
	public static String generateRandomSequence(String[] nucleotides, float[] nucleotProb, int n){
		StringBuilder shuffledSeq= new StringBuilder();
		int i= 1;
		// Random rnd= new Random();
		while(i <= n){
			float r = rnd.nextFloat();
			float cdf = (float) 0.0;
			for(int z = 0; z < nucleotides.length; z++){
				cdf += nucleotProb[z];
				if (r <= cdf){
					shuffledSeq.append(nucleotides[z]);
					break;
				}
			}
			i++;
		}
		return(shuffledSeq.toString());
	}	
	
	/**
	 * Produce a Markov model as defined in package markovChain.
	 * @param order Order of the markov chain
	 * @param sequence Sequence to model
	 * @param pctUnobserved Assign pseudocounts up to this percentage to unobserved
	 * kmers.
	 * @return Markov model
	 */
	public static markovChain.Model generateMarkovModel(
			int order, String sequence, double pctUnobserved){
		
		markovChain.Alphabet alphabet= new markovChain.Alphabet();
		Set<Character> alphaset= alphabet.getAlphabetFromString(sequence);
		
		markovChain.Model model= new markovChain.Model(order, alphaset, sequence, pctUnobserved);
		return(model);
	}
	
	/**
	 * Generate a random sequence of length n using the input Markov model.
	 * @param model
	 * @param n
	 * @return
	 */
	public static String generateRandomSequence(markovChain.Model model, int n){
		markovChain.MarkovSequence mc= new markovChain.MarkovSequence(); 		
		String rndSequence= mc.generateMarkovSequence(model, n);
		return(rndSequence);
	}
	
	public static String reverseComplement(String sequence){
		DNASequence dnaseq= new DNASequence(sequence, ambiguityCompoundSet);
		String revcomp= dnaseq.getReverseComplement().getSequenceAsString();
		return( revcomp );
	}

	public static List<SeqMatcher> regexMatcher(Pattern pattern, String sequence, boolean revcomp){
		List<SeqMatcher> matchList= new ArrayList<SeqMatcher>();
		
		Matcher match= pattern.matcher(sequence);
	
		while (match.find()){
			matchList.add(new SeqMatcher(sequence, match, true));
		}
		if (revcomp){
			match= pattern.matcher(reverseComplement(sequence));
			while (match.find()){
				matchList.add(new SeqMatcher(sequence, match, false));
			}			
		}
		return(matchList);
	}

	public static int regexCounter(Pattern pattern, String sequence, boolean revcomp){
		int count= 0;
		Matcher match= pattern.matcher(sequence);
		while (match.find()){
			count++;
		}
		if (revcomp){
			match= pattern.matcher(reverseComplement(sequence));
			while (match.find()){
				count++;
			}			
		}
		return(count);
	}

	/**
	 * Map the observed number of hits to the null distribution to return the empirical pvalue
	 * for the obs. no. of hits coming from random chance.  
	 * @param nullDistr
	 * @param observedRegexHits
	 * @return
	 */
	public static float mapObservedHitsToNull(TreeMultiset<Integer> nullDistr, int observedRegexHits){
		// TreeSet<Integer> sortedHits= new TreeSet<Integer>(nullDistr.elementSet());
		int rankOfObsHits= 0;
		for(Integer nhits : nullDistr.elementSet()){
			int count= nullDistr.count(nhits);
			if (observedRegexHits > nhits){
				rankOfObsHits += count;
			}
			else {
				break;
			}
		}
		float empiricalPvalue= 1 - (rankOfObsHits / (float) nullDistr.size());
		return(empiricalPvalue);
	}
	
	public static Writer printer(String filename){
		
		Writer wr= null;
		if (filename == null || filename.equals("-")){
			wr= new OutputStreamWriter(System.out);
		}
		else {
			try {
				wr = new PrintWriter(filename);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}
		return(wr);
	}
	
	/**
	 * Format pvalue to string with number of decimals
	 * given by the number of digits in the number of 
	 * executed randomizations.
	 * @param p
	 * @param nrand
	 */
	public static String formatPvalue(float p, int nrand){
		p = Math.round(p * nrand) / (float)nrand;
		String fp = String.format("%." + String.valueOf(nrand).length() + "f", p);
		return(fp);
	}
	
	public static void pregressReport(long t0, long t1, int nwind, Window window, ProbsLookUpTable tab){
		StringBuilder sb= new StringBuilder();
		sb.append(window.getChrom()); sb.append(":"); sb.append(window.getStart()); sb.append("-"); sb.append(window.getEnd()); sb.append(";\t");
		sb.append(t1 - t0); sb.append(" ms\t");
		sb.append(nwind); sb.append(" windows\t");
		sb.append(tab.probsLookUpTable.size()); sb.append(" size look up;");
		String log= sb.toString();
		System.err.println(log);
	}
} // End class

