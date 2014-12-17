package sequenceMatcher;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.zip.GZIPInputStream;

public class Opener {
	
	/**
	 * Method to open a file given its name. 
	 * @param filename Name of file to open. Can be gzipped (must end in .gz) 
	 * if "-" then read from stdin (System.in).
	 * @return
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public static InputStream open(String filename) throws IOException{
		
		InputStream inStream;
		if(filename.equals("-")){
			inStream= System.in;
		} else if(filename.endsWith(".gz")){
			inStream = new GZIPInputStream(new FileInputStream(filename));
		} else {
			inStream = new FileInputStream( filename );
		}

		return inStream;
	}
	
	/**
	 * Open a file, possibly gzipped, or from stdin (recognized by "-") and return a 
	 * BufferedReader ready to iterate through.
	 * @param filename
	 * @return
	 * @throws IOException
	 */
	public static BufferedReader openBr(String filename) throws IOException{

		BufferedReader br;
		
		if(filename.equals("-")){
			br = new BufferedReader(new InputStreamReader(System.in));
		} else if(filename.endsWith(".gz")){
			InputStream gzipStream = new GZIPInputStream(new FileInputStream(filename));
			Reader decoder = new InputStreamReader(gzipStream);
			br= new BufferedReader(decoder);
		} else {
			br= new BufferedReader(new FileReader(filename));
		}
		
		return br;
		
	}

	
}
