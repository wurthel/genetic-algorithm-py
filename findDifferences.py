template = """
AGYDLLGDGRPETLWLGIGTLLMLIGTFYFLVRGWGVTDKDAREYYAVTILVPGIASAAYLSMFFGIGLTEVTVGGEMLDIYYARYADWLFTTPLLLLDLALLAKVDRVTIGTLVGVDALMIVTGLIGALSHTAIARYSWWLFSTICMIVVLYFLATSLRSAAKERGPEVASTFNTLTALVLVLWTAYPILWIIGTEGAGVVGLGIETLLFMVLDVTAKVGFGFILLRSRAILGDTEAPE
"""
sequences = """
AGYDLLGDGRPETLWLGIHTLLMSIGTFYFLVRGWGVTDKDAREYYAVTILQGGHNSAAFLSIFFGIGLTEVTVGGEMLDIHYARYADWLFTTPLLHLTLALLAKVDRVTIGTLVGVSAWMIVTGLWGSLEHTAIARYHWWLFSTICNIVNLYFLATSLRSAAKERGPEVASTFNTLEALVAVLWTALWILWWEGTLGAGVVGLGIETLLFMDLDVTTKVGFGFILLRSRAILGDTEAPE	601.8
AGYDLLGDGRPETLWLGIGTLLMLIGTFYFLVRGWGVTDKDAREYYAETILVPGDASAAYLHMFFGIGLTEVTVGGEMLSIYYTRYADWLFTIKLLRLHLALLAKVDRVTIGTLVEVDIHMIVTGTIGALSHTAIAMYSQWMFSTQCMEVHLYFLATSLRSAAKERGPEVASTFNTLTALVPVLWTAEPILWIIKTTGAGVVGLGIDTLLFMVLDVDAKSGFGFIKLRSRAILGDTEAPE	479.3
AGYDLLGDGRPETLWLGIGTLLMLIGTFYELVRGWGVTDKDAREYYAEPILVPRIASAAYLSMFFGIGLTEVTVGGEMLSIYYARYAEVLFTLPVLLLDLALLAKVDRVTIGTLVGVDALMIVTDLKGSLSHTAIARYSSWLFSTICLIKVLYFLATSLRSAAKERGPEVASTFNTLTADVWVLWTAYGQLWTIGTEGRGVVGLGIWTLLFMVLSVTAKVGFGFELLRSRAILGDTEAPE	436.1
AGYDLLGDGRPETLWLGRPTLLCLIDTFYFLVRGWGVTDKDAREYYAVTRLSPGIASAAYLKFFFGIGLTEVTVGGEMLDIYYARYADWLFTTPLLGLDLALLAKVDRVTIGTLVGVDALMIVTGLIGALSHTRIVRYSWWLFVHIQSIVVLYFLATSLRSAAKERGPEVASTFNTLTARTWVDWIAYPIQWYIGTQGVGVVGLGIPTNLFFVLDFTSKVAEGFILLRSRAILGDTEAPE	543.3
AGYDLLGDGRPETLWLGEYGLLMEIGIFYELVRGWGVTDKDAREYYARTRLVPGIQQDAILSWFFGIGLTEVTVGGEMLIYYYARYANWLFVTPWLLLDLALLAKVDRVTIGTLCGVDALMIVTGVIGALDHTAIQRYSHWTFSTNCMRVVLYFLATSLRSAAKERGPEVASTFNTLTALVAVLWTSYPETHIIGTEGAGVVGLGIETKLEMVSDVTCKVGFGFIKLRSRAILGDTEAPE	472.6
AGYDLLGDGRPETLWLGIGTLLMLIKTFYFLVRGWGVTDKDAREYYARAILVPGIWSAADLQMFFGIGLTEVTVGGEMLDIHYAYYDAYLHTTPLLLLKLALLAKVDRVTIGTLVGVDALTIVTDLAWALEHTHIARYSLQTFSTICWIPVLYFLATSLRSAAKERGPEVASTFNTLHALVLVLWTASPKQWIIHTEGEGVVGLGIETKLFMEADVTAKVRCGEILLRSRAILGDTEAPE	604.3
AGYDLLGDGRPETLWLGDGCLLMLIGTFYFLVRGWGVTDKDAREYYAVEFLCPGRTSKAYRSMFFGIGLTEVTVGGEMLDIRYVRYADWLFTTPELLLLLALLAKVDRVTIGTRVGKVALTCVTGLIGALSHTAIARYSWHLFSTICMIRVLYFLATSLRSAAKERGPEVASTFNTLTDLVLVLWTIMHILYIEGTEGQGVVGLGINDLLFTVLAVRAKVLEGFILLRSRAILGDTEAPE	549.4
AGYDLLGDGRPETLWLGIGELLVLIGTFYFLVRGWGVTDKDAREYYAVTILVPGKYSAYYDSMFFGIGLTEVTVGGEMLDIYYARLSDWEFTTGLYFRCLALLAKVDRVTIGTLVEVDALYIVTGWIGALSHTAIARYSWWLFFTICMIVHLYFLATSLRSAAKERGPEVASTFNTLTAKYFDLQTAAQILWIIRTEGAGVVGLGIETMQFMVLDVTAKVKFGFILLRSRAILGDTEAPE	515.1
AGYDLLGDGRPETLWLGIGTLLMLIGTFYFLVRGWGVTDKDAREYYAVTILVPGIASAAYLSMFFGIGLTEVTVGGEMLDIYYARYADWLFTTPLLLLDLALLAKVDRVTIGTLVGVDALMIVTGLIGALSHTAIARYSWWLFSTICMIVVLYFLATSLRSAAKERGPEVASTFNTLTALVLVLWTAYPILWIIGTEGAGVVGLGIETLLFMVLDVTAKVGFGFILLRSRAILGDTEAPE	562.5
AGYDLLGDGRPETLWLGHYTLLMLITTFYRLVRGWGVTDKDAREYYAVTILVPGIASFAYLSMFFGIGLTEVTVGGEMLQHMYAEYADCHFTTHLLLLDLALLAKVDRVTIGTLYGVDALMIVTGLIGALSHTEIAPKSVWLFSTIPMIVVLYFEATSLRSAAKERGPEVASTFNTLTALYLVLWTAYPILWDQGTIGAGVVGLGIETVLFMVLEVTAKVPSGDILLRSRAILGDTEAPE	540.3
AGYDLLGDGRPETLWLGIGTLLMLIKTFYFLVRGWGVTDKDAREYYARAILVPGIWSHADLYMFFGIGLTEVTVGGEMLDIHYAYNDAYLHTFPLLLLKLALLAKVDRVTIGTLVGYDALTCVTDYAWALEHTHIARYSFQTFSEICWIPVLYFLATSLRSAAKERGPEVASTFNTAHNHVLVLWEASPKQWSIHTEGEGVVGLGIETKLFPEADVTAKVRCGEILLRSRAILGDTEAPE	628.7
AGYDLLGDGRPETLWLGRPMLLVLIGTFYFLVRGWGVTDKDAREYYAVDIDVIGKYSDYYDSMFFGIGLTEVTVGGEMLDIYYARLSDWDFTTGLYFRCLALLAKVDRVTIGTNVEVDALNEVTTWIGALSHTAIARYSWHLAFTTCMIGHLYFKATSLRSAAKERGPEVASTFNTLQAKYFDWQTAANILWCIWTEGAGVVGLGIYTEQFMVLGVTAKKWFGKIKLRSRAILGDTEAPE	507.7
AGYDLLGDGRPETLWLGIGKLLCLIGTFYFLVRGWGVTDKDAREYYAVTILVPGSASAAYLDMFFGIGLTEVTVGGEMLDIYYARYAHWLPTLELLLLKLALLAKVDRVTIGTLVGVDALMIVTGLEGALFHTAIARYSWWLKVTICWYRYLYFLATSLRSAAKERGPEVASTFNTDTALVLVLWWAYPILWMIGTEGAGVVGLGIQTLRCHVMDVIAKEGKGFDLLRSRAILGDTEAPE	436.2
AGYDLLGDGRPETLWLGIPTLLDSIGTFYYLVRGWGVTDKDAREYYAVTDLQGGHNVIEYLSDFFGIGLTEVTVGGEMLDIHYVWYADWLFLTKLLICTRALLAKVDRVTIGTEVHVSNEMTVTGLWGSLWHTHIARYHYWLRQTICNIVNLYFLATSLRSAAKERGPEVASTFNTLEALVAFLWTACWILWWNGTLGAGVVGLGIETLLSMDLDVITKHGDGEIDLRSRAILGDTEAPE	453.5
AGYDLLGDGRPETLWLGHYTLLMLIQLFYRLVRGWGVTDKDAREYYADTILEPGDADFAYTSMFFGIGLTEVTVGGEMLQRYYAKYADCHGTGHLLLLILALLAKVDRVTIGTLYRVNHLMIVTTLQGALSHTIIAPRSGWWFYPHPMHVELYFEATSLRSAAKERGPEVASTFNTLTARYLGSWTDYWILTDQGTIGAGVVGLGIDTVLFMVDMVTAKVPSGDIRLRSRAILGDTEAPE	441.5
AGYDLLGDGRPETLWLGIGELLVLIGTFYFLVRGWGVTDKDAREYYAVTRLVPGKYSAYYDSMFFGIGLTEVTVGGEMLDIYYARLSDWEFTAGLYFCCLALLAKVDRVTIGTLVEVTALNIVTGSIGANRHTAIARYSWWDFATFCMIVHLYFNATSLRSAAKERGPEVASTFNTLTAKYFDLQTAAQIWWSIRTLGAGVVGLGIPTMQFMVLEVTQKVKFGFIHLRSRAILGDTEAPE	578.7
AGYDLLGDGRPETLWLGIGTLLMLIKTFYFLVRGWGVTDKDAREYYARAILVPGIWPHADLYVFFGIGLTEVTVGGEMLDIHYDYNDAYAHTTPLLLLKWALLAKVDRVTIGTLVGYDALTCVTDYAWQLKHTHIALDSAQTESEICWIPDLYFLATSLRSAAKERGPEVASTFNTAENHVLVLWEASPKQWNIHTEGWGVVGLGIHTTLFEGADVTAKVVCGEILLRSRAILGDTEAPE	584.8
AGYDLLGDGRPETLWLGKYGLLMEIGIFYTLVRGWGVTDKDAREYYAQTRLVPSCQQDAIDGWFFGIGLTEVTVGGEMLDYKYDSMAPWLFVTPWLVLDLALLAKVDRVTIGTLCTVHAMMIVTGQIGALDHTSIQRYSHWTCSTHCMRVVLYFEATSLRSAAKERGPEVASTFNTLTALVAVLWTSYPETLIISTEGAGVVGLGIEYKLMMKSDVHCKLGFGFIKLRSRAILGDTEAPE	436.1
AGYDLLGDGRPETLWLGWGCLLMRIGTFYFLVRGWGVTDKDAREYYAVMFDCPLRTSKAYRSMFFGIGLTEVTVGGEMLYIYYVHYADWLFTTPELLLLDALLAKVDRVTIGTRVQQVALTCVTGLHQAHSHTAIADYHWELFSTICHIRVLYFLATSLRSAAKERGPEVASTFNTDTDLVLVLWADSHILYIEGTEGQGVVGLGISDLLFTVLDVYNKHLEGFIELRSRAILGDTEAPE	520.8
AGYDLLGDGRPETLWLGIMELLVLIGEFYFLVRGWGVTDKDAREYYAVTRLVPGKYSAYYDSMFFGIGLTEVTVGGEMLYIEYACLKMQEFTAELYFRQLALLAKVDRVTIGTKVLVTVLNIVVQTIQTNRHTAIARYSWWMRATFCMDVHLYFNATSLRSAAKERGPEVASTFNTMTCKYFDLQTAGQISWSIRTLGAGVVGLGIPTCRFHVLEVTQKVHFGFIHLRSRAILGDTEAPE	511.1
AGYDLLGDGRPETLWLGYHGLLMEIRIFYRLVRGWGVTDKDAREYYAQWRLVPSCQQYANDGWFFGIGLTEVTVGGEMLPYKYKQMAPWLFVTPHLKLKLALLAKVDRVTIGTLCEVTWRMIVTYIIGALDHTSIALYEPWTCSTHCMRVVLYFEATSLRSAAKERGPEVASTFNTLTADVDVLWTSYPFTLIIFTYGAGVVGLGIEYKLYTKSDVHCKLGFGFISLRSRAILGDTEAPE	537.1
AGYDLLGDGRPETLWLGIGTLKMLIGSFYDLVRGWGVTDKDAREYYACAILVKGKWPHADLYRFFGIGLTEVTVGGEMLWIHYPYNIAYAHTMPLLRLKWALLAKVDRVTIGTLVGYDALTCVTAYAWYKSHTDIALDSAPTTAEICWIPDLYFEATSLRSAAKERGPEVASTFNTAENHVLKLWEASPKYWTIHTCAHGVVGLGIITTLTEGADVTAKVVCGAILLRSRAILGDTEAPE	609.9
AGYDLLGDGRPETLWLGIMELKSLIGEFYFLVRGWGVTDKDAREYYAVTRLVPGKYSMYYDSMFFGIGLTEVTVGGEMLYHLYDCLKMQEFDAELMFAQSALLAKVDRVTIGTKVLVNVLNIVSQTDQVNRHTQIACYSWWMRATFCMAVHLYFLATSLRSAAKERGPEVASTFNTMTSDYFDLQTAGQHSWLIRTRGAGVVGLGIPTCRFHVSETTMKVHFGFIHLRSRAILGDTEAPE	466.7
AGYDLLGDGRPETLWLGKFKLLMEIGIFYGLVRGWGVTDKDAREYYAQTRLVDSCQQDAIDGDFFGIGLTEVTVGGEMLDYTYDQMAPCLFVTPWLRLCHALLAKVDRVTIGTLCWYHAMMIVTGQDGALDHTSIQRYSPWTCPTHCCRVVLYFEATSLRSAAKERGPEVASTFNTLTALVAVLWTSYQETLIEITTGAGVVGLGIEYKLMMKSDVHIKLGFGFRWLRSRAILGDTEAPE	482.3
AGYDLLGDGRPETLWLGIMELLHLIGEFYVLVRGWGVTDKDAREYYAVTRLYPYKYSAYYDSVFFGIGLTEVTVGGEMLYIGYACLKMQEGAAELFFRQLALLAKVDRVTIGTKVLDTVLNIVVQMIQTNRHTAIAGYSWWMRATFCMDVHLYFNATSLRSAAKERGPEVASTFNTMTCKVFDLQHAGYISWSIRTLGCGVVGLGIPTCYFHVLEVRQKVYFGFIHLRSRAILGDTEAPE	530.4
AGYDLLGDGRPETLWLGMGCLLMLIGTFYELVRGWGVTDKDAREYYAVEFLFPGRPSKAMRSHFFGIGLTEVTVGGEMLDFHYVEYADWLFTTPYTLLLDALLAKVDRVTIGTHVGKVALTCVLGLISALSHTKIAEYSWHLFSTIPMIRVLYFLATSLRSAAKERGPEVASTFNTLTDVVLVLWTIMHILYIEGTHGHGVVGLGICKLLFKVTAVRAKVLEGFKLLRSRAILGDTEAPE	591.1
AGYDLLGDGRPETLWLGCGELLVLIGTFYFLVRGWGVTDKDAREYYAKTRLVKGKYSAAVDSEFFGIGLTEVTVGGEMLDIYYACLMPDNVTAGLQFCNLALLAKVDRVTIGTLVEVTALNTVTGCDGADGHTAIARVSWWDFATFEMIWHLYFNATSLRSAAKERGPEVASTFNTLTAKYFDLQTAAQIWWSIRTFIEGVVGLGIPTMQFMVLEVTQKVKFGFIHLRSRAILGDTEAPE	668.8
"""

template = template.strip()
sequences = sequences.split()
results = []
for i in range(0, len(sequences), 2):
    sequence, v = sequences[i].strip(), float(sequences[i+1])
    oline = ""
    N = 0
    for j in range(len(template)):
        if template[j] != sequence[j]:
            oline += f"{template[j]}/{j+1}/{sequence[j]}, "
            N += 1
    results.append((oline, v, f"Count {N}"))

results = sorted(results, key=lambda x: x[1], reverse=True)
fout = open("differences.txt", "w")
for x in results:
    print(x[0], " ", x[1], " ", x[2], file=fout)
fout.close()