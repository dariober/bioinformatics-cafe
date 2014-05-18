package markovChain;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Sets;

public class Model {
	
	private Map<String, HashMultiset<Character>> actualModel;
	private Map<String, HashMultiset<Character>> emptyModel;
	private Map<String, HashMultiset<Character>> adjustedModel;
	private int order;
	private Set<Character> alphabet;
//	private int actualModelSize;
	private double pctUnobserved; 	// The percentage assigned to the unobserved 
									// combinations of kmers and letters.
									// Add pseudocounts to the unobserved elements
									// to avoid the MC getting stuck.

	/* --------------------------------------------------
	 *        S E T T E R S  &  G E T T E R S
	 * ----------------------------------------------- */
	
	public Map<String, HashMultiset<Character>> getActualModel() {
		return actualModel;
	}

	public void setActualModel(Map<String, HashMultiset<Character>> actualModel) {
		this.actualModel = actualModel;
	}

	public Map<String, HashMultiset<Character>> getAdjustedModel() {
		return adjustedModel;
	}

	public void setAdjustedModel(Map<String, HashMultiset<Character>> adjustedModel) {
		this.adjustedModel = adjustedModel;
	}
	
	public int getOrder() {
		return order;
	}

	public void setOrder(int order) {
		if (order < 1){
			System.err.println("Order must be >= 1. Got: " + order);
			System.exit(1);
		}
		this.order = order;
	}
	
	public Set<Character> getAlphabet() {
		return alphabet;
	}

	public void setAlphabet(Set<Character> alphabet) {
		this.alphabet = alphabet;
	}
	
	public Map<String, HashMultiset<Character>> getEmptyModel() {
		return emptyModel;
	}

	public void setEmptyModel(Map<String, HashMultiset<Character>> emptyModel) {
		this.emptyModel = emptyModel;
	}

	public double getPctUnobserved() {
		return pctUnobserved;
	}

	public void setPctUnobserved(double pctUnobserved) {
		if (pctUnobserved < 0 || pctUnobserved >= 1){
			System.err.println("Percentage unobserved elements 0 =< x < 1. Got: " + pctUnobserved);
			System.exit(1);			
		}
		this.pctUnobserved = pctUnobserved;
	}
	
	/* --------------------------------------------------
	 *              C O N S T R U C T O R S
	 * ----------------------------------------------- */
	
	public Model(){		
	} 
	
	public Model(int order, Set<Character> alphabet, String sequence, double pctUnobserved){
		if (order < 1){
			System.err.println("Order must be >= 1. Got: " + order);
			System.exit(1);
		}
		if (pctUnobserved < 0 || pctUnobserved >= 1){
			System.err.println("Percentage unobserved elements 0 =< x < 1. Got: " + pctUnobserved);
			System.exit(1);			
		}
		this.order= order;
		this.pctUnobserved= pctUnobserved;
		this.alphabet= alphabet;
		this.generateEmptyModel();
		this.getModelFromSequence(sequence);
		this.adjustModelForMissingObs();
	}
	
	/**
	 * Generate model for sequence given a preexisting model.
	 * Useful to avoid re-creating empty models every time the same
	 * order and alphabet is going to be used. 
	 * @param model
	 * @param sequence
	 */
	public Model(Model model, String sequence){
		this.order= model.order;
		this.pctUnobserved= model.pctUnobserved;
		this.alphabet= model.alphabet;
		this.generateEmptyModel();
		this.getModelFromSequence(sequence);
		this.adjustModelForMissingObs();
	}
	
	/* --------------------------------------------------
	 *                    M E T H O D S
	 * ----------------------------------------------- */
	
	/**
	 * Initialize a makov model by generating all the possible kmers. Each kmer
	 * assigned a count of 1 of each nucleotide.
	 * @param nucleotides Set of nucleotides to build kmers. Typically A, C, G, T
	 * @param order Order of the Markov chain. I.e. Length of each kmer
	 * @return Map of kmer and initialized counts.
	 */
	private void generateEmptyModel(){
	
		ImmutableList.Builder<Set<Character>> iList= ImmutableList.builder();
		int i= 0;
		while(i < order){
			iList.add(alphabet);
			i++;
		}
		Set<List<Character>> kmerSet= Sets.cartesianProduct(iList.build());
	
		Map<String, HashMultiset<Character>> model= new HashMap<String, HashMultiset<Character>>();
		for(List<Character> kmerAsList : kmerSet){
			// List of chars to string as key for map
			StringBuilder sb= new StringBuilder();
			for (char c : kmerAsList){
				sb.append(c);
			}
			String kmer= sb.toString();
			HashMultiset<Character> nucCount= HashMultiset.create();
			for (char x : alphabet){
				// Initilaize to 1 each nucleotide
				nucCount.add(x);
			}
			model.put(kmer, nucCount);
		}
		this.emptyModel= model;
	}

	/**
	 * Edit the input sequence to cope with characters not in the model
	 * alphabet.
	 * @param sequence Sequence to handle
	 * @param mode What to do with unexpected characters. Allowed modes are
	 * 	- "pass": Do nothing, return the string as is.
	 *  - "rm": Remove unexpected chars.
	 *  - "replace": Replace them with the character given in repl.
	 * @param repl: Replace unexpected chars with this. Used or not depending on mode. 
	 * @return Sequence edited accordingly
	 */
	public String handleUnexpSeqChars(String sequence, String mode, char repl){
		if (mode.equals("pass")){
			return(sequence);
		}
		else if (mode.equals("rm")){
			StringBuilder sb= new StringBuilder();
			for (char c : sequence.toCharArray()){
				if (this.alphabet.contains(c)){
					sb.append(c);
				}
			}
			if (sb.length() < this.order){
				System.err.println("Length of sequence after removing unexpected chars "
						+ "(" + sb.length() + ") is less then the order of the Markov chain (" + order + ")");
				System.exit(1);
			}
			return(sb.toString());
		}
		else if (mode.equals("replace")){
			StringBuilder sb= new StringBuilder();
			for (char c : sequence.toCharArray()){
				if (this.alphabet.contains(c)){
					sb.append(c);
				}
				else {
					sb.append(repl);
				}
			}
			return(sb.toString());
		}
		else {
			System.err.println("Unexpected keyword for arg mode: " + mode);
			System.exit(1);
			return("Wrong");
		}
	}
	
	/**
	 * Create a model by scanning input string to record and count all the 
	 * kmers of length = "order".
	 * @param sequence
	 * TODO: Handle characters not in the alphabet.
	 */
	private void getModelFromSequence(String sequence){
		sequence= handleUnexpSeqChars(sequence, "rm", ' ');

//		sequence= sequence + '\0';
		
		Map<String, HashMultiset<Character>> model= 
				new HashMap<String, HashMultiset<Character>>();
		
		for (int i= 0; i < (sequence.length() - order); i++){
			String kmer= sequence.substring(i, i + order);
			// System.out.println(kmer);
			char nextLetter= sequence.charAt(i + order);
			if (!model.containsKey(kmer)){
				HashMultiset<Character> nucCount= HashMultiset.create();
				nucCount.add(nextLetter);
				model.put(kmer, nucCount);
			} 
			else {
				HashMultiset<Character> nucCount= model.get(kmer);
				nucCount.add(nextLetter);
				model.put(kmer, nucCount);
			}
		}
		this.actualModel= model;
	}
	
	/**
	 * Set the number of items in the models before adjustment for 
	 * missing elements.
	 * This number be should equal to the length of the sequence used to
	 * build the model minus the order (length of kmer).
	 */
	public int modelSize(Map<String, HashMultiset<Character>> model){
		int modelSize= 0;
		for (String k : model.keySet()){
			modelSize += model.get(k).size();
		}
		return(modelSize);
	}
	
	public int countModelElements(Map<String, HashMultiset<Character>> model){
		int cnt= 0;
		for (String kmer : model.keySet()){
			cnt += model.get(kmer).elementSet().size();
		}
		return(cnt);
	}
	
	/**
	 * 
	 * MEMO: The counts to redistribute come from this equation:
	 *
	 *                              countMissingElements
	 *	pctUnobserved =  --------------------------------------------
     *            		 countMissingElements + sum(actualCounts) + X
	 *
	 * 	X= countsToRedistribute
	 *  The countsToRedistribute are partioned among the actual counts 
	 *  proportinally to the actual count weight.
	 *
	 * R CODE TO REPRODUCE THIS FUNCTION:
	 * ----------------------------------  
	    pctUnobserved<- 0.3
		actualModelSize<- 300
		actualCounts<- rep(20, 15)
		fullModelSize<- 16      
		
		countMissingElements<- fullModelSize - length(actualCounts)
		countsToRedistribute<- (countMissingElements - (pctUnobserved * (countMissingElements + actualModelSize))) / pctUnobserved
		if (countsToRedistribute < 0){
		    countsToRedistribute<- 0
		}
		
		stopifnot(
		    countMissingElements / (countsToRedistribute + actualModelSize + countMissingElements) <= pctUnobserved
		)
		
		countWeights<- (actualCounts / sum(actualCounts))
		actualCountsAdjusted<- actualCounts + (countsToRedistribute * countWeights)
		
		adjustedAndPseudocounts<- c(actualCountsAdjusted, rep(1, countMissingElements))
		
		stopifnot(
		    countMissingElements / sum(adjustedAndPseudocounts) <= pctUnobserved
		)
	 */
	private void adjustModelForMissingObs(){
		if (pctUnobserved == 0){
			// Effectively disable model adjustment
			this.adjustedModel= actualModel;
			return;
		}
		double actualElementsCount= countModelElements(this.actualModel);
		double actualModelSize= modelSize(this.actualModel);
		double fullModelSize= Math.pow((double)alphabet.size(), (double)order) * (double)alphabet.size(); 

		double countMissingElements= fullModelSize - actualElementsCount;
		double countsToRedistribute= (countMissingElements - (pctUnobserved * (countMissingElements + actualModelSize))) 
				/ pctUnobserved;
		if (countsToRedistribute < 0){
			// Reset to zero since we only augment the actual counts or leave as they are, never decrease.
			countsToRedistribute= 0;
		}
		double ratio= countMissingElements / (countsToRedistribute + actualModelSize + countMissingElements);
		if (ratio * 0.99 > pctUnobserved){
			// The ratio countMissingElements / total counts after adjustment must be <=  the required
			// pctUnobserved. Allow for some rounding error.
			System.err.println("Invalid computation of adjusted counts");
			System.err.println("ratio: " + ratio);
			System.err.println("actualElementsCount: " + actualElementsCount);
			System.err.println("actualModelSize: " + actualModelSize);
			System.err.println("fullModelSize: " + fullModelSize);
			System.err.println("countMissingElements: " + countMissingElements);
			System.err.println("countsToRedistribute: " + countsToRedistribute);
			System.exit(1);
		}
		
		// FOR DEBUGGING:
		// System.out.println("actualElementsCount:" + actualElementsCount);
		// System.out.println("actualModelSize:" + actualModelSize);
		// System.out.println("fullModelSize:" + fullModelSize);
		// System.out.println("countMissingElements:" + countMissingElements);
		// System.out.println("countsToRedistribute:" + countsToRedistribute);
		
		// Augment the actual counts proportionally to their weight
		// to redistribute the pseudocounts in countsToRedistribute
		Map<String, HashMultiset<Character>> adjustedModel= 
				new HashMap<String, HashMultiset<Character>>();
		
		for (String kmer : this.actualModel.keySet()){
			HashMultiset<Character> actualCounts= this.actualModel.get(kmer); 
			
			// Create a multiset with all the letters (nucleotides).
			HashMultiset<Character> adjustedCounts= HashMultiset.create();
			adjustedCounts.addAll(alphabet);
			
			// Now put the adjusted counts
			for (char c : actualCounts.elementSet()){
				double actualCount= actualCounts.count(c);
				double countWeight= actualCount / actualModelSize;
				double actualCountAdjusted= actualCount + (countsToRedistribute * countWeight);
				if (actualCountAdjusted > 2147483647){
					System.err.println("Adjusted counts greater then int size: " + actualCountAdjusted);
				}
				adjustedCounts.setCount(c, (int) Math.round(actualCountAdjusted));
			}
			adjustedModel.put(kmer, adjustedCounts);
		}
		// Fill up the adjustedCounts with the unobserved elements
		for (String kmer : emptyModel.keySet()){
			if (!adjustedModel.containsKey(kmer)){
				HashMultiset<Character> pseudoNucCounts= emptyModel.get(kmer);
				adjustedModel.put(kmer, pseudoNucCounts);
			}
		}
		if ((countMissingElements / modelSize(adjustedModel) * 0.9) > pctUnobserved){
			System.err.println("Missing elements make up more than the required percentage: " 
					+ countMissingElements / modelSize(adjustedModel));
			System.err.println("actualModel: " + actualModel);
			System.err.println("adjustedModel: " + adjustedModel);
			System.err.println("adjustedModelSize: " + modelSize(adjustedModel) );
			System.exit(1);
		}
		this.adjustedModel= adjustedModel;
	}
	
	public String toString(){
		
		StringBuilder sb= new StringBuilder();
		sb.append("MC Order: ");		sb.append(order); 			sb.append("\n");
		sb.append("Alphabet: ");		sb.append(alphabet); 		sb.append("\n");
		sb.append("Actual model: ");	sb.append(actualModel);		sb.append("\n");
		sb.append("Elements in actual model: ");	sb.append(countModelElements(actualModel));		sb.append("\n");
		sb.append("Empty model: ");		sb.append(emptyModel);		sb.append("\n");
		sb.append("Adjusted model: ");	sb.append(adjustedModel);	sb.append("\n");
		sb.append("Actual model size: ");		sb.append(modelSize(this.actualModel));	sb.append("\n");
		sb.append("Adjusted model size: ");		sb.append(modelSize(this.adjustedModel));
		return(sb.toString());
	}

}
