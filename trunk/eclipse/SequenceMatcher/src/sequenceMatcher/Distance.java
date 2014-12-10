package sequenceMatcher;

public class Distance {

	/**
	 * Modified from http://stackoverflow.com/questions/16260752/using-for-loop-to-get-the-hamming-distance-between-2-strings
	 * @param sequence1
	 * @param sequence2
	 * @param nmax Maximum number of edits allowed. Returns -1 if 
	 * 	exceeded. Set to -1 to allow any number of edits.
	 * @return
	 */
	public static int hammingDist(String sequence1, String sequence2, int nmax) {
	    
		// Corner cases
		boolean same= sequence1.equals(sequence2);
		if(same){
			return 0;
		}
		if(!same && nmax == 0){
			return -1;
		}
		// -----------------------------
		char[] s1 = sequence1.toCharArray();
	    char[] s2 = sequence2.toCharArray();

	    int shorter = Math.min(s1.length, s2.length);
	    int longest = Math.max(s1.length, s2.length);

	    int result = 0;
	    for (int i=0; i<shorter; i++) {
	        if (s1[i] != s2[i]) {
	        	result++;
	        	if(nmax >= 0 && result > nmax){
	        		break;
	        	}
	        }
	    }

	    result += longest - shorter;
	    if(nmax >= 0 && result > nmax){
	    	result= -1;
	    }
	    return result;
	}
		
}

/**
 * Compute Levenshtein distance between string s0 and s1
 * Code taken from 
 * http://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#Java
 * @param s0
 * @param s1
 * @return
 */
/*
public static int levenDist(String s0, String s1) {                          
    int len0 = s0.length() + 1;                                                     
    int len1 = s1.length() + 1;                                                     
 
    // the array of distances                                                       
    int[] cost = new int[len0];                                                     
    int[] newcost = new int[len0];                                                  
 
    // initial cost of skipping prefix in String s0                                 
    for (int i = 0; i < len0; i++) cost[i] = i;                                     
 
    // dynamicaly computing the array of distances                                  
 
    // transformation cost for each letter in s1                                    
    for (int j = 1; j < len1; j++) {                                                
        // initial cost of skipping prefix in String s1                             
        newcost[0] = j;                                                             
 
        // transformation cost for each letter in s0                                
        for(int i = 1; i < len0; i++) {                                             
            // matching current letters in both strings                             
            int match = (s0.charAt(i - 1) == s1.charAt(j - 1)) ? 0 : 1;             
 
            // computing cost for each transformation                               
            int cost_replace = cost[i - 1] + match;                                 
            int cost_insert  = cost[i] + 1;                                         
            int cost_delete  = newcost[i - 1] + 1;                                  
 
            // keep minimum cost                                                    
            newcost[i] = Math.min(Math.min(cost_insert, cost_delete), cost_replace);
        }                                                                           
 
        // swap cost/newcost arrays                                                 
        int[] swap = cost; cost = newcost; newcost = swap;                          
    }                                                                               
 
    // the distance is the cost for transforming all letters in both strings        
    return cost[len0 - 1];                                                          
}
*/