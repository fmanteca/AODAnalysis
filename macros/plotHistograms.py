""" -> PLOT THE HISTOGRAMS <-

>> Command: 

   python plotHistograms.py --inputFile output.root --outputDir plots_output --log 

"""


import ROOT
import argparse



#################################################################################

def tuneHistogram(histo, color):

    histo.SetFillColorAlpha(color, 0.3)
    histo.SetLineColor(color)
    histo.GetXaxis().SetLabelSize(0.03)
    histo.GetYaxis().SetLabelSize(0.03)
    histo.GetXaxis().SetTitleSize(0.037)
    histo.GetYaxis().SetTitleSize(0.037)
    histo.GetXaxis().SetTitleOffset(1.27)
    histo.GetYaxis().SetTitleOffset(1.4)
 

#################################################################################
#################################################################################
#################################################################################


if __name__ == '__main__':

    ########## main code here
    parser = argparse.ArgumentParser(description = "Receive the parameters")
    parser.add_argument('--inputFile', action = 'store', type = str, dest = 'inputFile', help = 'Input ROOT file with the data')
    parser.add_argument('--outputDir', action = 'store', type = str, dest = 'outputDir', default = '', help = 'Name of the output dir')
    parser.add_argument('--log', action = 'store_true', dest = 'log', help = 'Linear orllogaritmic scale in Y axis')
    parser.add_argument('--norm', action = 'store', type = str, dest = 'norm', default = 'norm0', help = 'Normalization Option')
    args = parser.parse_args()


    ########## Open the tree

    file_input = ROOT.TFile(args.inputFile)

    samples = []
    variables = []

    ## List of samples and variables:
    for key in file_input.GetListOfKeys():
    
        h = key.ReadObj()
        histo_name = str(h.GetName())
        variables.append(histo_name) 


    ########## Load the colors
    colors = {}
    exec(open('colors.py', 'r'))
    
    for key in colors.keys():
        samples.append(key)


    ########## Plot the histograms

    c1 = ROOT.TCanvas('c1', '', 500, 500)
    c1.SetTicks(1,1)
    c2 = ROOT.TCanvas('c2', '', 500, 500)
    c2.SetTicks(1,1)
    c3 = ROOT.TCanvas('c3', '', 500, 500)
    c3.SetTicks(1,1)

    ROOT.gStyle.SetOptStat(0)
    ROOT.TGaxis.SetMaxDigits(4)
    prefix = '' # default

    for var in variables:

        ## Y axis scale options

        if (args.log): 
            prefix = 'log_'
            c1.SetLogy(1)
        else:
            prefix = 'lin_'
            c1.SetLogy(0)


        ## Legend construction

        leg = ROOT.TLegend(0.54, 0.89 - 3*0.04, 0.89, 0.89)
        leg.SetTextSize(0.03)
        leg.SetBorderSize(0)


        ## Overlap representation implementation

        c = 0

        for sam in samples:
            
            #h = file_input.Get(sam + '_' + var)
            h = eval('file_input.' + var)
            tuneHistogram(h, colors[sam])

            if (args.norm == 'norm1'):
                h.Scale(1/h.Integral())


            if (args.log): h.SetMaximum(10.*h.GetMaximum())
            else: h.SetMaximum(1.2*h.GetMaximum())

            if c == 0: 
                h.Draw('hist')
            else: 
                h.Draw('hist same')
                h.Draw('axis same')


            if not 'genpt' in var:
                leg.AddEntry(h, sam, 'f')



            c+=1
            
            if var.startswith('Muon1'):
                h2 = eval('file_input.' + var.replace('Muon1', 'Muon2'))
                h2.Add(h)
                tuneHistogram(h2, colors[sam])
                
                if 'genpt' in var:

                    if 'genpt' in var:
                        leg.AddEntry(h, "Gen_pt", 'f')

                    txt = var.replace('genpt', 'TuneP_pt')
                    h4 = eval('file_input.' + var.replace('genpt', 'TuneP_pt'))
                    h5 = eval('file_input.' + txt.replace('Muon1', 'Muon2'))
                    h5.Add(h4)
                    h5.SetLineColor(3)
                    h5.SetLineWidth(3)
                    leg.AddEntry(h5, 'TuneP_pt', 'f')
                    leg.Draw()

            else:
                c1.cd()
                h.Draw()
                leg.Draw()

                c1.SaveAs(str(args.outputDir) + prefix + str(var) + '.png')
                c1.SaveAs(str(args.outputDir )+ prefix + str(var) + '.pdf')
                
                c1.Clear()


        ## Draw the legend:
        




        ## Save the plots



        if var.startswith('Muon1'):

            c2.cd()
            if (args.log): 
                prefix = 'log_'
                c2.SetLogy(1)

            h2.Draw()
            leg.Draw()
            c2.SaveAs(str(args.outputDir) + prefix + str(var) + '_add_Muon2.png')
            c2.SaveAs(str(args.outputDir )+ prefix + str(var) + '_add_Muon2.pdf')
            
            c2.Clear()

            if 'genpt' in var:
                c3.cd()
                if (args.log): 
                    prefix = 'log_'
                    c3.SetLogy(1)
                    
                h2.SetMaximum(13000)
                h2.Draw()
                h5.Draw("SAME")
                leg.Draw()
                c3.SaveAs(str(args.outputDir) + prefix + str(var) + '_add_Muon2_same_TunePpt.png')
                c3.SaveAs(str(args.outputDir )+ prefix + str(var) + '_add_Muon2_sameTunePpt.pdf')
                
                c3.Clear()
                
                
            
