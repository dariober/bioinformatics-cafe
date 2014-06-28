package bisReadBias;

import java.io.*;


public class ScriptRunner {

	private int returnCode;
	private String stdout;
	private String stderr;

	/*              S E T T E R S     &      G E T T E R S                   */
	
	public String getStdout() {
		return stdout;
	}

	public void setStdout(String stdout) {
		this.stdout = stdout;
	}

	/**
	 * @return Stdout from Rscript
	 */
	public String getStderr() {
		return stderr;
	}

	public void setStderr(String stderr) {
		this.stderr = stderr;
	}

	/**
	 * @return Stderr from Rscript
	 */
	public int getReturnCode() {
		return returnCode;
	}

	public void setReturnCode(int returnCode) {
		this.returnCode = returnCode;
	}
	
	/*                           M E T H O D S                               */
	
	/**
	 * Read the resource file and write it to a temp file for later
	 * use.
	 * @param rscript
	 * @return Path to newly written temp file readable by others.
	 * @see http://stackoverflow.com/questions/20389255/java-reading-a-resource-file-from-within-jar
	 */
	public String getTmpFileFromResourceFile(String rscript){
		
		//Open input file as stream
		BufferedReader br= null;
		String outfile= null;
		try {
			// Read script file
			br = new BufferedReader(
					new InputStreamReader(getClass().getResourceAsStream(rscript)));
		
			// Prepare tmp output
			File temp = File.createTempFile("plotter_", ".tmp.R");
			temp.deleteOnExit();
			
			// Write script to tmp
			FileWriter fw= new FileWriter(temp);
			String line;
			while ((line= br.readLine()) != null){
				fw.write(line + "\n");
			}
			br.close();
			fw.close();
			
			// Return full path to tmp
			outfile= temp.getAbsolutePath();
			
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		return(outfile);
	}

	/**
	 * Run Rscript to execute "script" with the given argument string 
	 * @param script Script to run
	 * @param args String of args passed to `Rscript script`
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void runRscript(String script, String args ) throws IOException, InterruptedException {
	
		Runtime r = Runtime.getRuntime();
		Process p= null;
		if (script.equals("")){
			p = r.exec("Rscript " + "--version");
		} else {
			p = r.exec("Rscript " + script + " " + args);
		}
		p.waitFor();
				
		this.returnCode= p.exitValue();
		
		// Print both stdout and stderr from Rscript
		BufferedReader bout = new BufferedReader(new InputStreamReader(p.getInputStream()));
		String line = "";
		while ((line = bout.readLine()) != null) {
			this.stdout+=line + "\n";
		}
		bout.close();
		
		BufferedReader berr = new BufferedReader(new InputStreamReader(p.getErrorStream()));
		line = "";
		while ((line = berr.readLine()) != null) {
			this.stderr+=line + "\n";
		}
		berr.close();
		
	}

}
