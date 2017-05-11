import os, Params as params
main_directory = os.path.dirname(os.path.realpath(__file__))
STYLE_PATH = main_directory + '/Style.css'
SCRIPT_PATH = main_directory + '/Script.js'

class Phylip_format(object):
    def __init__(self, path, sample_names):
        self.path = path
        self.sample_names = sample_names
        
    def write_sites(self, sites):
        f = open(self.path,'w')
        f.write("cells\t")
        for site in sites:
            f.write(site.CHROM + ":" + str(site.real_POS()) + "\t")

        for sample in self.sample_names:
            f.write("\n") 
            f.write(sample + "\t")
            for site in sites:
                if (site.samples[sample].info == "HET" or site.samples[sample].info == "ADO-A1"):
                    f.write("1\t")
                elif (site.samples[sample].info == "HOMO-R" or site.samples[sample].info == "ADO-R"):
                    f.write("0\t")
                else:
                    f.write("?\t")
class HTML(object):
    def __init__(self, path, sample_names):
        self.path = path
        self.source_code = self.open()
        self.sample_names = sample_names

    def open(self):
        source_code = ['<html>\n<head>\n<style>\n']
        source_code.append(open(STYLE_PATH).read())
        source_code.append("</style>\n<script>\n")
        source_code.append(open(SCRIPT_PATH).read())
        source_code.append("</script>\n</head>\n<body>")
        source_code.append('<body><div id="button"> <button type="button" id="btn" onclick="change(\'min\')">Zoom out</button> <img height="30px" src="https://s17.postimg.org/bb0oek5i7/conbase.gif" style="float:right;  margin-right:20px"><br/>')
        params_line = 'dp_ms_limit: ' + str(params.dp_ms_limit) + ', msp_ratio: ' + str(params.msp_ratio) + ', msp_internal_ratio: ' + str(params.msp_internal_ratio) + '</div>' 
        source_code.append(params_line)
        source_code.append('\n<table id="data">\n')
        return source_code
    
    def write_site(self, site):
        self.source_code.append('<tr onclick="show(' + "'hidden_" +  site.CHROM + "_" + str(site.real_POS()) + "'" + ')">\n')

        self.source_code.append('<td class="info cell">chr' + site.CHROM + ':' + str(site.real_POS()) + '</td>\n')  
        self.source_code.append('<td class="bulk cell">' + site.get_bulk_info() + '</td>')  


        for sample_name in self.sample_names:
            self.source_code.append('<td class="'+ site.samples[sample_name].info  + ' cell" >' + sample_name + '\n')
            self.source_code.append('<table class="hidden ' + 'hidden_' +  site.CHROM + '_' + str(site.real_POS()) + '">\n')
            self.source_code.append('<tr><td colspan=5>' + site.samples[sample_name].get_AD(site) + '</td></tr>')
            self.source_code.append('<tr><td colspan=5>-----------</td></tr>\n')
            for snp_pos, msp in site.samples[sample_name].MSP.items():
                if msp.get_ms_total() > 0:
                    if msp.voted:
                        if site.snp_ms_win[snp_pos] is None:
                            self.source_code.append('<tr><td class="voted" colspan=5>SNP: ' + str(snp_pos+1) + '</td></tr>\n')
                        else:
                            self.source_code.append('<tr><td class="voted" colspan=5>SNP: ' + str(snp_pos+1) + ' (' + str(site.snp_ms_win[snp_pos]) + ')</td></tr>\n')
                    else:
                        self.source_code.append('<tr><td colspan=5>SNP: ' + str(snp_pos+1) + '</td></tr>\n')

                    self.source_code.append('<tr><td colspan=2>GT distribution</td>\n')

                    self.source_code.append('<td>' + str(msp.homo_R) + '</td>\n')
                    self.source_code.append('<td>' + str(msp.het) + '</td>\n')
                    self.source_code.append('<td>' + str(msp.homo_A1) + '</td></tr>\n')
                    self.source_code.append('<tr><td>MSP</td>\n')
                    self.source_code.append('<td> RR: ' + str(msp.RR) + '</td>\n')
                    self.source_code.append('<td> RA: ' + str(msp.RA) + '</td>\n')
                    self.source_code.append('<td> AR: ' + str(msp.AR) + '</td>\n')
                    self.source_code.append('<td> AA: ' + str(msp.AA) + '</td></tr>\n')
                    self.source_code.append('<tr><td colspan=5>-----------</td></tr>\n')
            self.source_code.append('</table>\n')
            self.source_code.append('</td>')
        self.source_code.append('</tr>\n')
       
        
    def close(self):
        self.source_code.append('\n</table>\n')
        self.source_code.append('</body>\n</html>')
        f = open(self.path,'w')
        f.write("".join(self.source_code))
        f.close()
        
class Site2JSON(object):
    def __init__(self, path):
        self.f = open(path,'w')    
    def write(self,site):
        if site.TYPE == '':
            self.f.write(repr(site) + '\n')

    def close(self):
        self.f.close()
