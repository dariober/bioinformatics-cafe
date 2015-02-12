package sequenceMatcher;


import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;
import org.biojava3.alignment.Alignments;
import org.biojava3.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava3.alignment.Alignments.PairwiseSequenceScorerType;
import org.biojava3.alignment.SimpleGapPenalty;
import org.biojava3.alignment.SubstitutionMatrixHelper;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;

import com.google.common.base.Joiner;

/**
 * Class to hold matches between two sequences.
 * nameA/B: Name of seqs
 * NM: NUmber of mismatches
 * strand: Whether match happens on revcomp(seqA) and seqB.
 * @author berald01
 *
 */
public class Match {

	/**
	 * HEADER is where you set what goes to output and in which order.
	 * The keys of Map are the positions of the columns in output.
	 * The values are used as header and to retrieve the appropriate
	 * attributes. See method toString(). 
	 * NB: Positions are 0-based.
	 */
	public static final Map<Integer, String> HEADER;
    static {
        LinkedHashMap<Integer, String> header= new LinkedHashMap<Integer, String>();
        header.put(0, "seq_A");
        header.put(1, "seq_B");
        header.put(2, "strand");
        header.put(3, "LD");
        header.put(4, "HD");
        header.put(5, "JWD");
        header.put(6, "len_A");
        header.put(7, "len_B");
        header.put(8, "pos");
        header.put(9, "NM");
        header.put(10, "aln_score");
        header.put(11, "n_ident");
        header.put(12, "n_sim");
        header.put(13, "pct_ident");
        header.put(14, "len_aln");
        header.put(15, "alnA");
        header.put(16, "alnB");
        HEADER= Collections.unmodifiableMap(header);
    }

    
	/**
	 * Reverse the HEADER map: Now keys are column names and values are the 
	 * positions.
	 */
	private static final Map<String, Integer> rev_header;
    static {
        LinkedHashMap<String, Integer> revheader= new LinkedHashMap<String, Integer>();
        for(int pos : HEADER.keySet()){
        	revheader.put(HEADER.get(pos), pos);
        }
        rev_header= Collections.unmodifiableMap(revheader);
    }

    
	private String nameA= ".";
	private String nameB= ".";
	private String strand= ".";
	private String seqA= ".";
	private String seqB= ".";
	private SequencePair<DNASequence, NucleotideCompound> alnPair= null;
	private Integer LD= -1; // Levenshtein dist
	private Integer HD= -1; // Hamming dist
	private Double JWD= -1.0; // JaroWinkler dist
	private int len_A;
	private int len_B;
	private int pos= -1; // Aln start pos on reference (seqA): Leftmost position as per sam spec
	private int NM= -1; // Nucleotide difference
	private int aln_score= -1; // Alignment score
	private int len_aln= -1;
	private int n_ident= -1;
	private int n_sim= -1;
	private double pct_ident= -1;
	private String alnA= ".";
	private String alnB= ".";
	private String alnMethod;

	/*         S e t t e r s    a n d   G e t t e r s  */
	public String getNameA() {
		return nameA;
	}
	public void setNameA(String nameA) {
		this.nameA = nameA;
	}
	public String getNameB() {
		return nameB;
	}
	public void setNameB(String nameB) {
		this.nameB = nameB;
	}
	public String getStrand() {
		return strand;
	}
	public void setStrand(String strand) {
		this.strand = strand;
	}
	public String getSeqA() {
		return seqA;
	}
	public void setSeqA(String seqA) {
		this.seqA = seqA;
	}
	public String getSeqB() {
		return seqB;
	}
	public void setSeqB(String seqB) {
		this.seqB = seqB;
	}

	public int getAln_score() {
		return aln_score;
	}
	
	public SequencePair<DNASequence, NucleotideCompound> getAlnPair() {
		return alnPair;
	}

	public String getAlnA() {
		return alnA;
	}
	public void setAlnA(String alnA) {
		this.alnA = alnA;
	}
	public String getAlnB() {
		return alnB;
	}
	public void setAlnB(String alnB) {
		this.alnB = alnB;
	}

	public double getPct_ident() {
		return pct_ident;
	}
	public void setPct_ident(double pct_ident) {
		this.pct_ident = pct_ident;
	}

	public int getPos() {
		return pos;
	}
	
	public int getNM() {
		return NM;
	}
	public void setNM(int nM) {
		NM = nM;
	}
	
	public int getLen_A() {
		return len_A;
	}
	public void setLen_A(int len_A) {
		this.len_A = len_A;
	}
	public int getLen_B() {
		return len_B;
	}
	public void setLen_B(int len_B) {
		this.len_B = len_B;
	}
	public void setLD(Integer lD) {
		LD = lD;
	}
	public int getLD() {
		return this.LD;
	}

	public void setHD(Integer hD) {
		HD = hD;
	}
	public int getHD() {
		return this.HD;
	}

	public void setJWD(Double jWD) {
		JWD = jWD;
	}
	public double getJWD(){
		return this.JWD;
	}

	public int getLen_aln() {
		return len_aln;
	}
	public void setLen_aln(int len_aln) {
		this.len_aln = len_aln;
	}

	public int getN_ident() {
		return n_ident;
	}
	public void setN_ident(int n_ident) {
		this.n_ident = n_ident;
	}
	public int getN_sim() {
		return n_sim;
	}
	public void setN_sim(int n_sim) {
		this.n_sim = n_sim;
	}
	public void setAlnPair(SequencePair<DNASequence, NucleotideCompound> alnPair) {
		this.alnPair = alnPair;
	}
	
	public String getAlnMethod() {
		return alnMethod;
	}
	public void setAlnMethod(String alnMethod) {
		this.alnMethod = alnMethod;
	}
	
	/* ------------------------------------------------------------------------ */
	// Distances
	// ---------
	public void computeJWD(){
		this.JWD= StringUtils.getJaroWinklerDistance(seqA, seqB);
	}	

	public void computeLD(int nmax) {
		if(nmax < 0){
			this.LD= StringUtils.getLevenshteinDistance(seqA, seqB);
		} else{
			this.LD= StringUtils.getLevenshteinDistance(seqA, seqB, nmax);
		}
	}
	public void computeLD() {
		this.LD= StringUtils.getLevenshteinDistance(seqA, seqB);
	}
	
	public void computeHD(int nmax) {
		this.HD= Distance.hammingDist(seqA, seqB, nmax);
	}
	public void computeHD() {
		this.HD= Distance.hammingDist(seqA, seqB, -1);
	}
	
	// -------------------------------------------------------------------------
	// C O N S T R U C T O R S
	Match(){
		
	}
	
	/**
	 * Initialize Match obj with two sequences. Sequences
	 * are represented by String arrays of length 2: Sequence name, Sequence. 
	 * @param seqA
	 * @param seqB
	 */
	Match(String[] seqA, String[] seqB){
		this.nameA= seqA[0];
		this.nameB= seqB[0];
		this.seqA= seqA[1];
		this.seqB= seqB[1];
		this.len_A= this.seqA.length();
		this.len_B= this.seqB.length();
	}
	
	// -------------------------------------------------------------------------
		
	/**
	 * Match sequence seqA with sequence seqB to get a distance
	 * metric using a specified method. This metric is going
	 * to be used for deciding whether additional metrics have to be computed
	 * and to decide whether the match is to be printed out.  
	 * @param method Method to compute mismatches
	 * @param nmax Max number of mismatches to output the match
	 * @return
	 */
	public double getFilterDistance(String method, int nmax){
		
		double d= -1;
		if(method.equals("LD")){
			this.computeLD(nmax);
			d= this.getLD();
		} else if (method.equals("HD")){
			this.computeHD(nmax);
			d= this.getHD();
		} else {
			System.err.println("Unsupported method: " + method);
			System.exit(1);
		}
		return d;
	}
	
	public void align(){
		
		AmbiguityDNACompoundSet ambiguityDNACompoundSet= new AmbiguityDNACompoundSet();
		
		SubstitutionMatrix<NucleotideCompound> matrix = SubstitutionMatrixHelper.getNuc4_4();	 
		SimpleGapPenalty gapP = new SimpleGapPenalty();
				
		
		PairwiseSequenceScorerType scorerType= null;
		PairwiseSequenceAlignerType alignerType= null;
		if(this.getAlnMethod() == null ){
			System.err.println("Alignment method not set!");
			System.exit(1);
		}
		if(this.getAlnMethod().equals("global")){
			scorerType= PairwiseSequenceScorerType.GLOBAL;
			alignerType= PairwiseSequenceAlignerType.GLOBAL;
		} else if (this.getAlnMethod().equals("local")) {
			scorerType= PairwiseSequenceScorerType.LOCAL;
			alignerType= PairwiseSequenceAlignerType.LOCAL;
		} else {
			System.err.println("Invalid option for alignment method: " + this.getAlnMethod());
			System.exit(1);
		}
		
		// Nomenclature Ours -> SAM -> Biojava:
		// seqA -> ref -> query
		// seqB -> read -> target
		DNASequence refQuery= new DNASequence(this.seqA, ambiguityDNACompoundSet);
		DNASequence readTarget= new DNASequence(this.seqB, ambiguityDNACompoundSet);
		
		int aln_score= Alignments.getAllPairsScores(
				Arrays.asList(refQuery, readTarget), scorerType, new SimpleGapPenalty(), matrix)[0];		
		
		SequencePair<DNASequence, NucleotideCompound> alnPair= Alignments.getPairwiseAlignment(
				refQuery, readTarget, alignerType, gapP, matrix);
		
		this.pos= alnPair.getIndexInQueryAt(1);		
		this.alnPair= alnPair;
		this.len_aln= alnPair.getLength();
		this.aln_score= aln_score;
		this.n_ident= alnPair.getNumIdenticals();
		this.n_sim= alnPair.getNumSimilars();
		this.pct_ident= (double) alnPair.getNumIdenticals() / (double) alnPair.getLength();
		this.alnA= this.alnPair.getAlignedSequences().get(0).getSequenceAsString();
		this.alnB= this.alnPair.getAlignedSequences().get(1).getSequenceAsString();
		this.NM= Distance.getNMtagFromAln(this.alnB, this.alnA);
	}	
	
	/*
	 * Print a match in tabular format. To change the order of the output change
	 * the keys in HEADER map. To add or remove change both here and HEADER.
	 * */
	public String toString(){

		String[] sa= new String[HEADER.size()];
		
		for(int k : HEADER.keySet()){
			String h= HEADER.get(k);
			if(h.equals("seq_A")) 		{sa[k]= this.nameA;} else
			if(h.equals("seq_B")) 		{sa[k]= this.nameB;} else
			if(h.equals("strand")) 		{sa[k]= this.strand;} else
			if(h.equals("LD")) 			{sa[k]= Integer.toString(this.LD);} else
			if(h.equals("HD")) 			{sa[k]= Integer.toString(this.HD);} else
			if(h.equals("JWD")) 		{sa[k]= String.format("%.2f", this.JWD);} else
			if(h.equals("len_A")) 		{sa[k]= Integer.toString(this.len_A);} else
			if(h.equals("len_B")) 		{sa[k]= Integer.toString(this.len_B);} else
			if(h.equals("pos")) 		{sa[k]= Integer.toString(this.pos);} else
			if(h.equals("NM")) 			{sa[k]= Integer.toString(this.NM);} else
			if(h.equals("aln_score")) 	{sa[k]= Integer.toString(this.aln_score);} else
			if(h.equals("n_ident"))   	{sa[k]= Integer.toString(this.n_ident);} else
			if(h.equals("n_sim")) 		{sa[k]= Integer.toString(this.n_sim);} else
			if(h.equals("pct_ident")) 	{sa[k]= String.format("%.2f", this.pct_ident);} else
			if(h.equals("len_aln")) 	{sa[k]= Integer.toString(this.len_aln);} else
			if(h.equals("alnA")) 		{sa[k]= this.alnA;} else
			if(h.equals("alnB")) 		{sa[k]= this.alnB;} 
			else {
				System.err.println("Unexpected header");
				System.exit(1);
			}
		}		
		return(Joiner.on("\t").join(sa));
	} 

	public String getHeader() {
		
		String[] sa= new String[HEADER.size()];
		for(int k : HEADER.keySet()){
			sa[k]= HEADER.get(k);
		}			
		return(Joiner.on("\t").join(sa));
	}
	
	/**
	 * Takes a string representing a line from the output of 
	 * SequenceMatcher match in tab format and returns a Match
	 * object from it.
	 * @param line
	 * @return
	 */
	public void stringToMatchObj(String line) {

		String[] aline= line.trim().split("\t", -1);
		if(aline.length != rev_header.size()){
			System.err.println("Incorrectly formatted string: " + aline.length + " filed found. Expected: " + rev_header.size());
			System.err.println("String: " + line);
			System.err.println("Header: " + this.getHeader());
			System.exit(1);
		}
		this.nameA= aline[rev_header.get("seq_A")];
		this.nameB= aline[rev_header.get("seq_B")];
		this.strand= aline[rev_header.get("strand")];
		this.LD= Integer.parseInt(aline[rev_header.get("LD")]);
		this.HD= Integer.parseInt(aline[rev_header.get("HD")]);
		this.JWD= Double.parseDouble(aline[rev_header.get("JWD")]);
		this.len_A= Integer.parseInt(aline[rev_header.get("len_A")]);
		this.len_B= Integer.parseInt(aline[rev_header.get("len_B")]);
		this.pos= Integer.parseInt(aline[rev_header.get("pos")]);
		this.NM= Integer.parseInt(aline[rev_header.get("NM")]);
		this.aln_score= Integer.parseInt(aline[rev_header.get("aln_score")]);
		this.n_ident= Integer.parseInt(aline[rev_header.get("n_ident")]);
		this.n_sim= Integer.parseInt(aline[rev_header.get("n_sim")]);
		this.pct_ident= Double.parseDouble(aline[rev_header.get("pct_ident")]);
		this.len_aln= Integer.parseInt(aline[rev_header.get("len_aln")]);
		this.alnA= aline[rev_header.get("alnA")];
		this.alnB= aline[rev_header.get("alnB")];
	}

}
