package aligner;

import java.io.*;
import java.util.zip.GZIPInputStream;

/**
 * No frills method to read fastq files optionally gzipped. Sequence and quality are 
 * on one line (not checked!). There is no sanity check for whether the fastq file is 
 * correctly formatted.
 * 
 * Basic usage:
 * BufferedReader br= FastqReader.openFastq("test_data/fk018.fq.gz");
 * String[] fqread= FastqReader.getNextRead(br);
 * 
 * @author berald01
 *
 */
public class FastqReader {
	
	public static BufferedReader openFastq(String fastq) throws FileNotFoundException, IOException{
		
		BufferedReader br;
		if(fastq.endsWith(".gz")){
		    br=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fastq))));
		} else {
		   br=new BufferedReader(new FileReader(fastq));
		}
		return(br);
	}

	/**
	 * @param br
	 * @return String array of length 4 with 
	 * 1) read name (w/o leading @)
	 * 2) Sequence
	 * 3) Comment (w/o leading +)
	 * 4) Quality
	 * @throws IOException
	 */
	public static String[] getNextRead(BufferedReader br) throws IOException {
		String[] fqread= new String[4];
		
		// Name
		String name= br.readLine();
		if (name == null){
			return(null);
		}
		fqread[0]= name.substring(1, name.length());
		// Sequence
		fqread[1]= br.readLine();
		// Comment
		String comment= br.readLine();
		fqread[2]= comment.substring(1, comment.length());
		// Quality
		fqread[3]= br.readLine();	
		return(fqread);
	}
	
	/*
	long nLine=-1L;
	String line;
	while((line=r.readLine())!=null){
	    nLine++;
	    if(nLine%4!=0 || !nameset.contains(line)) continue;
	    System.out.println(line);
	    System.out.println(r.readLine());
	    System.out.println(r.readLine());
	    System.out.println(r.readLine());
	    nLine+=3;
	}
	r.close();
	*/
}
