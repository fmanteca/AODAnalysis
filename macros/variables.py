""" File to especify the variables that are going to be plotted 

>> The structure should appear as follows:
We have a variable dict: variable = {}
variable['variable1_key'] = {'name' : 'name in the file or expression',
                             'range' : [0, 1],
                             'label' : 'label'}



>> NOTE: Expression syntax
Before each variable it is required to write event. such as
event.PV_vx - event.RefittedPV_vx


"""

variables['nMuons'] = {'name' : 'event.nMuons',
                       'range' : [0, 4],
                       'nbins' : 5,
                       'label' : 'nMuons',
                       'SafeCuts': '1'
}

variables['Muon_TuneP_pt1'] = {'name' : 'event.Muon_TunePTrack_pt[0]',
                            'range' : [200, 3500],
                            'nbins' : 100,
                            'label' : 'Muon_TuneP_pt1 [GeV]',
                            'SafeCuts' : 'event.nMuons > 0',
}

variables['Muon_TuneP_pt2'] = {'name' : 'event.Muon_TunePTrack_pt[1]',
                            'range' : [200, 3500],
                            'nbins' : 100,
                            'label' : 'Muon_TuneP_pt2 [GeV]',
                            'SafeCuts' : 'event.nMuons > 1',
}

variables['Muon_TuneP_eta1'] = {'name' : 'event.Muon_TunePTrack_eta[0]',
                            'range' : [200, 3500],
                            'nbins' : 100,
                            'label' : 'Muon_TuneP_eta1',
                            'SafeCuts' : 'event.nMuons > 0',
}

variables['Muon_TuneP_eta2'] = {'name' : 'event.Muon_TunePTrack_eta[1]',
                            'range' : [200, 3500],
                            'nbins' : 100,
                            'label' : 'Muon_TuneP_eta2',
                            'SafeCuts' : 'event.nMuons > 1',
}

variables['Muon_genpt1'] = {'name' : 'event.Muon_Genpt[0]',
                            'range' : [200, 3500],
                            'nbins' : 100,
                            'label' : 'GenMuon_pt1 [GeV]',
                            'SafeCuts' : 'event.nMuons > 0',
}

variables['Muon_genpt2'] = {'name' : 'event.Muon_Genpt[1]',
                            'range' : [200, 3500],
                            'nbins' : 100,
                            'label' : 'GenMuon_pt2 [GeV]',
                            'SafeCuts' : 'event.nMuons > 1',
}

variables['Hit1_Muon1_GeomDet1_x'] = {'name' : 'event.Hit_x[0]',
                                      'range' : [-200, 200],
                                      'nbins' : 100,
                                      'label' : 'Hit_Muon1_GeomDet1_x [cm]',
                                      'SafeCuts' : '(event.nMuons > 0 and event.Muon_nGeomDets > 0 and event.Muon_nHits > 0)',
}

variables['Prop_Muon1_GeomDet1_x'] = {'name' : 'event.Prop_x[0]',
                                      'range' : [-200, 200],
                                      'nbins' : 100,
                                      'label' : 'Prop_Muon1_GeomDet1_x [cm]',
                                      'SafeCuts' : '(event.nMuons > 0 and event.Muon_nGeomDets > 0)',
}




