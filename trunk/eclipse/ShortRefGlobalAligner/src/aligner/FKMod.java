package aligner;

import java.util.Arrays;
import java.util.List;


public class FKMod {

	// Indexes are 1-based!!
	public static final List<Integer> INDEX_FIRST_BARCODE= Arrays.asList(new Integer[] {1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
	public static final List<Integer> INDEX_SECOND_BARCODE= Arrays.asList(new Integer[] {48, 49, 50, 51, 52, 53, 54, 55, 56, 57});
	public static final List<Integer> INDEX_MODIFICATIONS= Arrays.asList(new Integer[] {35, 43});
	
	public String getQueryBasesAtTargetIndexes(Aligner aln, List<Integer> targetIndexes) {

		StringBuilder queryBases= new StringBuilder();
		
		for (int alnIdx=1; alnIdx <= aln.getPair().getLength(); alnIdx++){
			
			int targetIdx= aln.getPair().getIndexInTargetAt(alnIdx);
		
			if(targetIndexes.contains(targetIdx)){
				
				String queryBase= aln.getPair().getCompoundInQueryAt(alnIdx).getBase();
				String targetBase= aln.getPair().getCompoundInTargetAt(alnIdx).getBase();
				if (!targetBase.equals("-")){
					queryBases.append(queryBase);
				}				
			}
		}
		
		if(queryBases.length() > targetIndexes.size()){
			System.err.println(queryBases);
			System.err.println(aln);
			System.exit(1);
		}
		return(queryBases.toString());
	}

}
