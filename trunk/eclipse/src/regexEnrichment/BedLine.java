package regexEnrichment;

public class BedLine {
	private String chrom;
	private int start;
	private int end;
	private String name;
	private String score; // pvalue for window enrichment
	private boolean strand= true;
	private int wstart; // Window start & end
	private int wend; 
	private String group; // Match group from regex
	
	BedLine(){
		
	}
	BedLine(Window w, SeqMatcher m, String s){
		this.chrom= w.getChrom();
		this.start= w.getStart() + m.getStart();
		this.end= w.getStart() + m.getEnd();
		this.name= w.getChrom() 
				+ "_" + w.getStart() + m.getStart() 
				+ "_" + w.getStart() + m.getEnd() 
				+ (m.isStrand() ? "_for" : "_rev");
		this.score= s;
		this.strand= m.isStrand();
		this.wstart= w.getStart();	
		this.wend= w.getEnd();
		this.group= m.getGroup();
	}
	
	public String toString(){
		String sep= "\t";
		StringBuilder sb= new StringBuilder();
		sb.append(chrom); sb.append(sep);
		sb.append(start); sb.append(sep);
		sb.append(end); sb.append(sep);
		sb.append(name); sb.append(sep);
		sb.append(score); sb.append(sep);
		sb.append((strand) ? "+" : "-"); sb.append(sep);
		sb.append(wstart); sb.append(sep);
		sb.append(wend); sb.append(sep);
		sb.append(group);
		return(sb.toString());
	}
}
