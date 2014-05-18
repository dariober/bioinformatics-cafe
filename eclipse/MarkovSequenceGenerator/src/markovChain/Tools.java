package markovChain;

import java.io.FileNotFoundException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.Writer;

public class Tools {
	
	public static Writer printer(String filename){	
		Writer wr= null;
		if (filename == null || filename.equals("-")){
			wr= new OutputStreamWriter(System.out);
		}
		else {
			try {
				wr = new PrintWriter(filename);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}
		return(wr);
	}
	
}
