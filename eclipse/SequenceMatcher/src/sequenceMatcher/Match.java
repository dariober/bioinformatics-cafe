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
	private Integer NM;
	private String strand= ".";
	private String seqA= ".";
	private String seqB= ".";
	private int NWscore; // Needleman-wunch alignment score
	private SequencePair<DNASequence, NucleotideCompound> alnPair;
	private Integer LD; // Levenshtein dist
	private Integer HD; // Hamming dist
	private double JWD; // JaroWinkler dist
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
	public int getNM() {
		return NM;
	}
	public void setNM(int nM) {
		NM = nM;
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

	// Distances
	
	public Integer getLD(int nmax) {
		int LD= StringUtils.getLevenshteinDistance(seqA, seqB, nmax);
		return LD;
	}
	public void setLD(Integer lD) {
		LD = lD;
	}
	public Integer getHD(int nmax) {
		int HD= Distance.hammingDist(seqA, seqB, nmax);
		return HD;
	}
	public void setHD(Integer hD) {
		HD = hD;
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
	
	public double getJaroWinklerDistance(){
		double jwd= StringUtils.getJaroWinklerDistance(seqA, seqB);
		return jwd;
	}
	
	/**
	 * Match sequence seqA with sequence seqB. 
	 * @param setA
	 * @param seqB
	 * @param method Method to compute mismatches
	 * @param nmax Max number of mismatches to output the match
	 * @param strand Match also to the strand.
	 * @return with elements
	 */
	public void getDistance(String method, int nmax){
		
		int d= -1;
		if(method.equals("Leven") && nmax >= 0){
			d= StringUtils.getLevenshteinDistance(seqA, seqB, nmax);
		} else if(method.equals("Leven") && nmax < 0){
			d= StringUtils.getLevenshteinDistance(seqA, seqB);
		} else if (method.equals("Hamming")){
			d= Distance.hammingDist(seqA, seqB, nmax);
		} else {
			System.err.println("Unsupported method: " + method);
			System.exit(1);
		}
		this.NM= d;
	}
	
	/**
	 * Get distance between strings using a valid method.
	 * @param method: Leven or Hamming
	 * @param nmax: Max number of mismatches. If exceeded -1 is returned
	 * @param revcomp: Should seqA be reverse complemented before matching?
	 */
	//public void getDistance(String method, int nmax, boolean revcomp){
	//	if(revcomp){
	//		AmbiguityDNACompoundSet ambDNAset= new AmbiguityDNACompoundSet();
	//		this.seqA= new DNASequence(this.seqA, ambDNAset).getReverseComplement().getSequenceAsString();
	//		this.getDistance(method, nmax);
	//		this.strand= "-";
	//	} else {
	//		this.getDistance(method, nmax);
	//		this.strand= "+";
	//	}
	//}
	
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
	
	}	
	
	public String alnToString(){

		String[] sa= new String[6];
		
		if(this.alnPair == null){
			for(int i= 0; i < sa.length; i++){
				sa[i]= ".";
			}
		} else {
			sa[0]= Integer.toString(this.NWscore);
			sa[1]= Integer.toString(this.alnPair.getLength());	  
			sa[2]= Integer.toString(this.alnPair.getNumIdenticals());
			sa[3]= Integer.toString(this.alnPair.getNumSimilars());
			sa[4]= this.alnPair.getAlignedSequences().get(0).getSequenceAsString();
			sa[5]= this.alnPair.getAlignedSequences().get(1).getSequenceAsString();
		}
		return(Joiner.on("\t").join(sa));

	} 

	public String toString(){
		
		String[] sa= new String[6];
		
		sa[0]= nameA; 
		sa[1]= nameB;
		sa[2]= strand;
		sa[3]= (this.NM == null) ? "." : (Integer.toString(NM));
		sa[4]= Integer.toString(seqA.length()); 
		sa[5]= Integer.toString(seqB.length());
		// sa[6]= seqA; 
		// sa[7]= seqB;
		
		return(Joiner.on("\t").join(sa) + "\t" + this.alnToString());
	} 

	public String getHeader() {
		String[] header= new String[12];

		header[0]= "seq_A"; 
		header[1]= "seq_B";
		header[2]= "strand";
		header[3]= "NM";
		header[4]= "len_A";
		header[5]= "len_B";
		header[6]= "aln_score";
		header[7]= "len_aln";
		header[8]= "n_ident";
		header[9]= "n_sim";
		header[10]= "aln_A";
		header[11]= "aln_B";
		return(Joiner.on("\t").join(header));
	}
	
	public void setHeader(String header) {
		this.header = header;
	}

	
}
