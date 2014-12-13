package sequenceMatcher;


import java.util.Arrays;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.biojava3.alignment.Alignments;
import org.biojava3.alignment.SimpleGapPenalty;
import org.biojava3.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava3.alignment.Alignments.PairwiseSequenceScorerType;
import org.biojava3.alignment.SubstitutionMatrixHelper;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.DNASequence;

import com.google.common.base.Joiner;

import org.biojava3.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;

/**
 * Class to hold matches between two sequences.
 * nameA/B: Name of seqs
 * NM: NUmber of mismatches
 * strand: Whether match happens on revcomp(seqA) and seqB.
 * @author berald01
 *
 */
public class Match {
	
	private String nameA= ".";
	private String nameB= ".";
	private String strand= ".";
	private String seqA= ".";
	private String seqB= ".";
	private int NWscore; // Needleman-wunch alignment score
	private SequencePair<DNASequence, NucleotideCompound> alnPair;
	private double pct_ident; 
	private String alnA;
	private String alnB;
	private Integer LD; // Levenshtein dist
	private Integer HD; // Hamming dist
	private Double JWD; // JaroWinkler dist
	private int NM; // Nucleotide difference
	@SuppressWarnings("unused")
	private String header;
	
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

	public int getNWscore() {
		return NWscore;
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

	public int getNM() {
		return NM;
	}
	public void setNM(int nM) {
		NM = nM;
	}
	
	/* ------------------------------------------------------------------------ */
	// Distances
	// ---------
	public void setJWD(){
		this.JWD= StringUtils.getJaroWinklerDistance(seqA, seqB);
	}	
	public double getJWD(){
		return this.JWD;
	}

	public void setLD(int nmax) {
		if(nmax < 0){
			this.LD= StringUtils.getLevenshteinDistance(seqA, seqB);
		} else{
			this.LD= StringUtils.getLevenshteinDistance(seqA, seqB, nmax);
		}
	}
	public void setLD() {
		this.LD= StringUtils.getLevenshteinDistance(seqA, seqB);
	}
	public int getLD() {
		return this.LD;
	}
	
	public void setHD(int nmax) {
		this.HD= Distance.hammingDist(seqA, seqB, nmax);
	}
	public void setHD() {
		this.HD= Distance.hammingDist(seqA, seqB, -1);
	}
	public int getHD() {
		return this.HD;
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
			this.setLD(nmax);
			d= this.getLD();
		} else if (method.equals("HD")){
			this.setHD(nmax);
			d= this.getHD();
		} else {
			System.err.println("Unsupported method: " + method);
			System.exit(1);
		}
		return d;
	}
	
	public void align(){
		
		AmbiguityDNACompoundSet ambiguityDNACompoundSet= new AmbiguityDNACompoundSet();
		
		DNASequence read= new DNASequence(); 
		read= new DNASequence(this.seqA, ambiguityDNACompoundSet);

		DNASequence reference= new DNASequence(this.seqB, ambiguityDNACompoundSet);
	
		SubstitutionMatrix<NucleotideCompound> matrix = SubstitutionMatrixHelper.getNuc4_4();	 
		SimpleGapPenalty gapP = new SimpleGapPenalty();
				
		List<DNASequence> lst= Arrays.asList(reference, read);
				
		int NWscore= Alignments.getAllPairsScores(lst,
				PairwiseSequenceScorerType.GLOBAL, new SimpleGapPenalty(), matrix)[0];
		
		SequencePair<DNASequence, NucleotideCompound> alnPair= Alignments.getPairwiseAlignment(
				read, 
				reference,
				PairwiseSequenceAlignerType.GLOBAL, gapP, matrix);
		
		this.alnPair= alnPair;
		this.NWscore= NWscore;
		this.pct_ident= (double) alnPair.getNumIdenticals() / (double) alnPair.getLength();
		this.alnA= this.alnPair.getAlignedSequences().get(0).getSequenceAsString();
		this.alnB= this.alnPair.getAlignedSequences().get(1).getSequenceAsString();
		this.NM= Distance.getNMtagFromAln(this.alnB, this.alnA);
	}	
	
	public String alnToString(){

		String[] sa= new String[7];
		
		if(this.alnPair == null){
			for(int i= 0; i < sa.length; i++){
				sa[i]= ".";
			}
		} else {
			sa[0]= Integer.toString(this.NWscore);
			sa[1]= Integer.toString(this.alnPair.getLength());	  
			sa[2]= Integer.toString(this.alnPair.getNumIdenticals());
			sa[3]= Integer.toString(this.alnPair.getNumSimilars());
			sa[4]= String.format("%.2f", this.pct_ident);
			sa[5]= this.alnA;
			sa[6]= this.alnB;
		}
		return(Joiner.on("\t").join(sa));

	} 

	public String toString(){
		
		String[] sa= new String[8];
		
		sa[0]= nameA; 
		sa[1]= nameB;
		sa[2]= strand;
		sa[3]= (this.LD == null) ? "." : (Integer.toString(LD));
		sa[4]= (this.HD == null) ? "." : (Integer.toString(HD));
		sa[5]= (this.JWD == null) ? "." : (Double.toString(JWD));
		sa[6]= Integer.toString(seqA.length()); 
		sa[7]= Integer.toString(seqB.length());
		
		return(Joiner.on("\t").join(sa) + "\t" + this.alnToString());
	} 

	public String getHeader() {
		String[] header= new String[15];

		header[0]= "seq_A"; 
		header[1]= "seq_B";
		header[2]= "strand";
		header[3]= "LD";
		header[4]= "HD";
		header[5]= "JWD";
		header[6]= "len_A";
		header[7]= "len_B";
		header[8]= "aln_score";
		header[9]= "len_aln";
		header[10]= "n_ident";
		header[11]= "n_sim";
		header[12]= "pct_ident";
		header[13]= "aln_A";
		header[14]= "aln_B";
		return(Joiner.on("\t").join(header));
	}
	
	public void setHeader(String header) {
		this.header = header;
	}

}
