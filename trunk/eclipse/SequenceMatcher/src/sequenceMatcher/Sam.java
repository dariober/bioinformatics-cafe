package sequenceMatcher;

import java.util.ArrayList;

import com.google.common.base.Joiner;

public class Sam {

	private String qname;
	private int flag;
	private String rname;
	private int pos;
	private int mapq= 255;
	private String cigar;
	private String rnext= "*";
	private int pnext= 0;
	private int tlen= 0;
	private String seq;
	private String qual= "*";
	private String tags;
	
	public String getQname() {
		return qname;
	}

	public void setQname(String qname) {
		this.qname = qname;
	}

	public int getFlag() {
		return flag;
	}

	public void setFlag(int flag) {
		this.flag = flag;
	}

	public void setFlag(String strand) {
		if(strand.equals("-")){
			this.flag = 16;
		} else {
			this.flag = 0;
		}
	}
	
	public String getRname() {
		return rname;
	}

	public void setRname(String rname) {
		this.rname = rname;
	}

	public int getPos() {
		return pos;
	}

	public void setPos(int pos) {
		this.pos = pos;
	}

	public int getMapq() {
		return mapq;
	}

	public void setMapq(int mapq) {
		this.mapq = mapq;
	}

	public String getCigar() {
		return cigar;
	}

	public void setCigar(String cigar) {
		this.cigar = cigar;
	}

	public String getRnext() {
		return rnext;
	}

	public void setRnext(String rnext) {
		this.rnext = rnext;
	}

	public int getPnext() {
		return pnext;
	}

	public void setPnext(int pnext) {
		this.pnext = pnext;
	}

	public int getTlen() {
		return tlen;
	}

	public void setTlen(int tlen) {
		this.tlen = tlen;
	}

	public String getSeq() {
		return seq;
	}

	public void setSeq(String seq) {
		this.seq = seq;
	}

	public String getQual() {
		return qual;
	}

	public void setQual(String qual) {
		this.qual = qual;
	}

	public String getTags() {
		return tags;
	}

	public void setTags(String tags) {
		this.tags = tags;
	}
	
	/* ------------------------------------------------------------------------ */
	
	public static String fastaListToSQHeader(ArrayList<String[]> fastaList) {
		StringBuilder sb= new StringBuilder();
		for(int i= 0; i < fastaList.size(); i++){
			sb.append("@SQ");
			sb.append("\tSN:");
			sb.append(fastaList.get(i)[0]);
			sb.append("\tLN:");
			sb.append(fastaList.get(i)[1].length());
			sb.append("\n");
		}
		return sb.toString().trim();
	}

	public static int getAlnStartPos(String read, String ref) {

		if(read.length() != ref.length()){
			System.err.println("Alignments of unequal length!");
			System.exit(1);
		}
		
		int start= 1;
		for(int i=0; i < read.length(); i++){
			char a= read.charAt(i);
			char r= ref.charAt(i);
			if(a != '-' && r != '-'){
				return start;
			} else if(r == '-'){
				continue;
			} else if(a == '-'){
				start++;
			} else {
				System.err.println("Invalid alignment or unexpected case!");
			}
		}
		return -1;
	}

	
	/**
	 * Return the type of cigar operation for the reference and read char.
	 * Possibilities:
	 * read: N -> M
	 * ref:  N
	 * 
	 * read: N -> I
	 * ref:  -
	 * 
	 * read: - -> D
	 * ref:  N
	 * 
	 * @param readChar
	 * @param refChar
	 * @return
	 */
	private static char getCigarOperation(char readChar, char refChar){
		char cigarOp= '\0';
		if(readChar != '-' && refChar != '-'){
			cigarOp= 'M';
		} else if(readChar != '-' && refChar == '-'){
			cigarOp= 'I';
		} else if(readChar == '-' && refChar != '-'){
			cigarOp= 'D';
		} else {
			System.err.println("Unexpected read and ref characters: " + readChar + " and " + refChar);
			System.exit(1);
		}
		return cigarOp;
	}

	/**
	 * Check the sum of cigar operator lengths equal the length of read. 
	 * SAM spec: "Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ"
	 * @param cigarOps cigar operator
	 * @param opLens operator lengths 
	 * @param alnRead Read as it has been passed to getCigarFromAln();
	 * @return
	 */
	private static boolean validateCigarLen(ArrayList<Character> cigarOps,
			ArrayList<Integer> opLens, String alnRead){
		// Get original read by stripping the insertions
		String read= getSeqFromAln(alnRead);
		int exLen= 0;
		for(int i= 0; i < cigarOps.size(); i++){
			char op= cigarOps.get(i);
			if(op == 'M' || op == 'I' || op == 'S' || op == '=' || op == 'X'){
				exLen += opLens.get(i);
			}
		}
		if(exLen == read.length()){
			return true;
		} else {
			return false;
		}
		
	}
	
	public static String getCigarFromAln(String read, String ref) {
		
		if(read.length() != ref.length()){
			System.err.println("Alignments of unequal length!");
			System.exit(1);
		}
		if(read.length() == 0){
			return "*";
		} 
		
		ArrayList<Character> cigarOps= new ArrayList<Character>();
		ArrayList<Integer> opLens= new ArrayList<Integer>();
		
		char prevOp= getCigarOperation(read.charAt(0), ref.charAt(0));
		int opLen= 1;
		char curOp;
		for(int i= 1; i < read.length(); i++){
			curOp= getCigarOperation(read.charAt(i), ref.charAt(i));
			if(curOp == prevOp){
				opLen++;
			} else {
				cigarOps.add(prevOp);
				opLens.add(opLen);
				prevOp= curOp;
				opLen= 1;
			}
		}
		cigarOps.add(prevOp);
		opLens.add(opLen);
		
		/*
		If the first or last operator is D, remove it as it doesn't make sense. E.g.
		--ACTG-- Read
		TTACTGCC Ref
		With 2D4M2D would mean you are claiming a deletion without having 
		sequenced it. 
		*/
		if(cigarOps.get(0) == 'D'){
			cigarOps.remove(0);
			opLens.remove(0);
		}

		if(cigarOps.get(cigarOps.size()-1) == 'D'){
			cigarOps.remove(cigarOps.size()-1);
			opLens.remove(opLens.size()-1);
		}

		if(cigarOps.size() != opLens.size()){
			System.err.println("Incorrect CIGAR string: Number of operators and lengths don;t match.");
			System.exit(1);
		}

		boolean ok= validateCigarLen(cigarOps, opLens, read);
		if(!ok){
			System.err.println("Sum of cigar operations does not equal read length.");
			System.exit(1);
		}
		
		StringBuilder sb= new StringBuilder();
		for(int i=0; i < cigarOps.size(); i++){
			sb.append(opLens.get(i));
			sb.append(cigarOps.get(i));
		}
		return sb.toString();
		
	}

	public static String getSeqFromAln(String alnRead) {
		String read= alnRead.replace("-", "");
		return read;
	}

	public static ArrayList<String> matchToTagList(Match m) {
		
		ArrayList<String> tagList= new ArrayList<String>();
		
		tagList.add("NM:i:" + Integer.toString(m.getNM()));
		tagList.add("AS:i:" + Integer.toString(m.getNWscore()));
		tagList.add("XL:i:" + Integer.toString(m.getLD()));
		tagList.add("XH:i:" + Integer.toString(m.getHD()));
		tagList.add("XJ:f:" + String.format("%.4f", m.getJWD()) );
		tagList.add("XP:f:" + String.format("%.4f", m.getPct_ident()) );		
		tagList.add("XR:Z:" + m.getAlnA() );
		tagList.add("XS:Z:" + m.getAlnB() );
		return tagList;
	}
	
	public static String matchToTagString(Match m){
		ArrayList<String> tagList= matchToTagList(m);
		String tagStr= Joiner.on("\t").join(tagList);
		return tagStr;
	}
	
	
	public static Sam matchToSam(Match m){
		// Memo: 
		// alnA -> "reference" -> RNAME 
		// alnB -> "read" -> QNAME
		
		Sam sam= new Sam();
		sam.setQname(m.getNameB());
		sam.setFlag(m.getStrand());
		sam.setRname(m.getNameA());
		sam.setPos(Sam.getAlnStartPos(m.getAlnB(), m.getAlnA()));
		sam.setMapq(255);
		sam.setCigar(Sam.getCigarFromAln(m.getAlnB(), m.getAlnA()));

		String seq= Sam.getSeqFromAln(m.getAlnB());
		String qual= new String(new char[seq.length()]).replace("\0", "I");
		sam.setSeq(seq);
		sam.setQual(qual);

		sam.setTags(Sam.matchToTagString(m));
		return sam;
	}
	
	public String toString(){
		// Here order matters!!
		StringBuilder sb= new StringBuilder();

		sb.append(qname); sb.append("\t"); //1
		sb.append(flag);  sb.append("\t"); //2
		sb.append(rname); sb.append("\t"); //3
		sb.append(pos);   sb.append("\t"); //4
		sb.append(mapq);  sb.append("\t"); //5
		sb.append(cigar); sb.append("\t"); //6
		sb.append(rnext); sb.append("\t"); //7
		sb.append(pnext); sb.append("\t"); //8
		sb.append(tlen);  sb.append("\t"); //9
		sb.append(seq);   sb.append("\t"); //10
		sb.append(qual);   sb.append("\t"); //11
		sb.append(tags);  sb.append("\t"); //12
		
		return sb.toString().trim();
	} 
}
