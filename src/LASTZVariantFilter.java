
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;


/**
 * Bespoke code for Andrew South.
 * 
 * Reads a tab delimited file of values and outputs a filtered, processed and extended version of that file.
 *
 * If the filepath is a directory then the output will be the concaternation of all *.interval files in that directory.
 * 
 * Notes from Andy:
 * 1)      The .interval file catalogs all variants detected from the sequence run for a particular sample. Because the sequence is error prone I am not interested in variants showing up only once, these are the majority. Therefore I opened the .variant file in excel and gave each column a header so that I could filter on the last column “No Variants”. I am not interested in variants that show up only once. In fact I only took 3 as the cut off but thinking now could you make this two please? This is the first worksheet of the excel file.
 * 2)      I then copied and pasted those rows with >3 (I’d like 2 for your script please) in the no variant column and made a new column “% reads”. Using this I filtered out those variants which, although were represented 3 or more times were less than 10% of the total reads. In the example I sent I am comparing “No Variants”, column K, with “No Reads”, column E, but I should have compared “No Qual”, column J instead of E, please could you do this. I would like the % printed as this gives me an allelic frequency which is useful.
 * 3)      I then copied and pasted those rows >10% into the next worksheet “filter” and this is where the manual nonsense came in. I had to go through and score the actual variants. The original variant output tells you what the base (A, C, G or T, columns G-J) should be, the reference base (normal if you like) but my next galaxy work flow needs the variant base. So I have had add a column “SNP” and then go through and put in the variant manually. As you can see this is usually one other base but sometimes there are two and so this has to be inputted A/C instead of A. Another way of doing this could be duplicating the row with different in each but your script ought to be able to generate A/C in such a circumstance. Additionally, I applied another filtering criteria here – if the variants were 2 bases instead of the normal 1 and then the split took each individual variant below the 10% I discarded them. For instance, 30 reads, 26 were G (reference) with 4 were variants (meeting initial criteria) but these variants were 2 A’s and 2 C’s, this takes each individual variant below the 10% cut-off so then I discard. This level of filtering is desired but not essential as I would check all detrimental variants after the next workflow anyway, so if this looks to be a pain forgo this particular filter.
 * 4)      From the filter worksheet in excel I then copied and pasted the data minus the headers into another worksheet “IC1A_SNP”, added another column at the end for the sample id and then saved this as a tab delimited txt file which is what I need the data in. I then put all these files together for all the samples into one large file. If your script could do a batch of these files, say in a directory and then join them up this would be such a time saving utility.
 * 
 * @author matt
 */
public class LASTZVariantFilter {
	static int FILTER_VARIANT_NUM=2;
	static double FILTER_PERCENT=10.0;
	static String FILE_EXTENSION=".interval";
	
	private String filepath = null;
	private int lines=0;
	private int processed=0;
	private List<Variant> list = new ArrayList<Variant>();
	
	/**
	 * This method is run automatically by Java if it exists.
	 * @param args see printHelp()
	 */
	public static void main(String[] args) throws FileNotFoundException, IOException {
		if (args.length>0) {
			if (args[0].equals("--help")) {
				printHelp();
			} else {
				String filepath = args[args.length-1];
				File file = new File(filepath);
				if (file.isDirectory()) {
					for (File child : file.listFiles()) {
						if (child.getName().endsWith(FILE_EXTENSION)) {
							runProcessor(child);
						}
					}
				} else {
					runProcessor(file);
				}
			}
		} else {
			printHelp();
		}
	}
	
	private static void runProcessor(File file) throws FileNotFoundException, IOException {
		LASTZVariantFilter processor = new LASTZVariantFilter();
		processor.setFilepath(file.getAbsolutePath());
		processor.summarise();
		for (Variant variant : processor.list) {
			if (variant.variants>=FILTER_VARIANT_NUM && variant.percentReads>=FILTER_PERCENT && variant.snp.length()>0) {
				System.out.println(variant.toString() + "\t" + file.getName().replace(FILE_EXTENSION, ""));
			}
		}												
	}

	private static void printHelp() {
		System.out.println("LASTZ Variant Filter");
		System.out.println("Requires java v6 or later (running \"java -version\" tells you the version).");
		System.out.println("\nUsage: java LASTZVariantFilter filepath");	
		System.out.println("\nIf filepath is a directory then all *.interval files within that directory will be");
		System.out.println("processed and the name of the originating file attached to the end of each row.");
	}
	
	/**
	 * Process the file
	 * @throws FileNotFoundException if the filepath provided in setFilepath() doesn't properly point to a file
	 * @throws IOException if something goes wrong while reading the file
	 */
	public void summarise() throws FileNotFoundException, IOException {
		File file = new File(filepath);
		FileReader reader = new FileReader(file);
		BufferedReader buff = new BufferedReader(reader);
		String line;
		while((line=buff.readLine())!=null) {
			if (line.length()>0) {
				lines++;
				String[] cols = line.split("\t");
				if (cols.length==11) {
					list.add(new Variant(cols));
					processed++;
				}
			}
		}
		buff.close();
		reader.close();
		// sense checking
		if (lines!=processed) {
			throw new RuntimeException("Error: Are you sure this is the right file? found " + lines + " lines but could only process " + processed);
		}
	}

	public void setFilepath(String filepath) {
		this.filepath=filepath;
	}

	public String getFilepath() {
		return filepath;
	}

	/** 
	 * Holds the values for a row of the input file and some extra calculations.
	 * 
	 * @author matt
	 *
	 */
	class Variant {
		String chr;
		long start;
		long end;
		String refBase;
		int reads;
		int a;
		int c;
		int g;
		int t;
		int qualReads;
		int variants;
		
		List<String> snpList=new ArrayList<String>();
		String snp;
		int maxVariants=-1;
		double percentReads;
		
		Variant(String[] line) {
			chr = line[0];
			start = Long.valueOf(line[1]);
			end = Long.valueOf(line[2]);
			refBase = line[3];
			reads = Integer.valueOf(line[4]);
			a = Integer.valueOf(line[5]);
			c = Integer.valueOf(line[6]);
			g = Integer.valueOf(line[7]);
			t = Integer.valueOf(line[8]);
			qualReads = Integer.valueOf(line[9]);
			variants = Integer.valueOf(line[10]);
			
			addSNPResult(a, "A");
			addSNPResult(c, "C");
			addSNPResult(g, "G");
			addSNPResult(t, "T");
			snp=snpList.toString().replace(", ", "/").replace("[","").replace("]","");
			if (snpList.size()>1) {
				percentReads=100.0*maxVariants/qualReads;
			} else {
				percentReads=100.0*variants/qualReads;				
			}
		}
		
		private void addSNPResult(int letter, String character) {
			if (letter>=FILTER_VARIANT_NUM && !refBase.toUpperCase().equals(character)) {
				snpList.add(character);
				if (letter>maxVariants) maxVariants=letter;
			}
		}
		
		public String toString() {
			return chr+"\t"+start+"\t"+end+"\t"+snp+"\t"+refBase+"\t"+a+"\t"+c+"\t"+g+"\t"+t+"\t"+qualReads+"\t"+variants+"\t"+percentReads;
		}
	}
}
