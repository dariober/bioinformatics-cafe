package markovChain;

import java.util.Map;
import java.util.Random;
import java.util.Set;

import com.google.common.collect.HashMultiset;

public class MarkovSequence {
	
	private static Random random= new Random();
	
	/**
	 * Sample a random kmer from a set of kmer's. Used to pick a random element to start 
	 * a Markov chain.
	 * Sampling weighted by the counts in the Multiset.
	 * @param model A model to pick a kmer from
	 * @return
	 */
	public String initializeChain(Map<String, HashMultiset<Character>> model){
		
		int modelSize= 0;
		for (String k : model.keySet()){
			modelSize += model.get(k).size();
		}
		int r= random.nextInt(modelSize);
		int cdf= 0;
		String outkmer= null;
		for (String kmer : model.keySet()){
			cdf += model.get(kmer).size();
			if (r < cdf){
				outkmer= kmer;
				break;
			}
		}
		return(outkmer);
	}
	

	/**
	 * Sample a character from markov model given a seed (kmer as key in the HashMultiset).
	 * See also http://www.realpython.com/blog/python/lyricize-a-flask-app-to-create-lyrics-using-markov-chains/#.UwfHrBhaX3g 
	 * @param model Markov model typically produced by generateModel.
	 * @param kmer String to query the model to get a character following this kmer.
	 * @return A character sampled from the options corresponding to this kmer. 
	 */
	public char getNextChar(Map<String, HashMultiset<Character>> model, String kmer){
		if (!model.containsKey(kmer)){
			System.err.println("\ngetNextChar: kmer '" + kmer + "' not found in model table.\n");
			System.exit(1);
		}
		HashMultiset<Character> nucCounts= model.get(kmer);
		//System.out.println(nucCounts);
		int rnd= random.nextInt(nucCounts.size());
		int cdf= 0;
		char nextChar= ' ';
		for (char x : nucCounts.elementSet()){
			cdf += nucCounts.count(x);
			if (rnd < cdf){
				nextChar= x;
				break;
			}
		}
		if (nextChar == ' '){
			System.err.println("\ngetNextChar: Could not sample a character.\n");
			System.exit(1);
		}		
		return(nextChar);
	}	
	
	public String generateMarkovSequence(Model model, int length){

		StringBuilder mcsequence= new StringBuilder();
		
		String kmer= initializeChain(model.getAdjustedModel());
//		System.out.print(kmer + " ");
		while(mcsequence.length() < length){
			char newChar= getNextChar(model.getAdjustedModel(), kmer);
//			System.out.print("kmer:" + kmer + " next:" + newChar + " ");
			mcsequence.append(newChar);
			kmer= kmer.substring(1, kmer.length()) + newChar;
		}
		return(mcsequence.toString());
	}	
	
}
