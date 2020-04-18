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

variables['Muon1_TuneP_pt'] = {'name' : 'event.Muon_TunePTrack_pt[0]',
                            'range' : [200, 3500],
                            'nbins' : 50,
                            'label' : 'Muon_TuneP_pt1 [GeV]',
                            'SafeCuts' : 'event.nMuons > 0',
}

variables['Muon2_TuneP_pt'] = {'name' : 'event.Muon_TunePTrack_pt[1]',
                            'range' : [200, 3500],
                            'nbins' : 50,
                            'label' : 'Muon_TuneP_pt2 [GeV]',
                            'SafeCuts' : 'event.nMuons > 1',
}

variables['Muon1_TuneP_eta'] = {'name' : 'event.Muon_TunePTrack_eta[0]',
                                'range' : [-2.4, 2.4],
                                'nbins' : 50,
                                'label' : 'Muon_TuneP_eta1',
                                'SafeCuts' : 'event.nMuons > 0',
                            }

variables['Muon2_TuneP_eta'] = {'name' : 'event.Muon_TunePTrack_eta[1]',
                                'range' : [-2.4, 2.4],
                                'nbins' : 50,
                                'label' : 'Muon_TuneP_eta2',
                                'SafeCuts' : 'event.nMuons > 1',
                            }

variables['Muon1_genpt'] = {'name' : 'event.Muon_Genpt[0]',
                            'range' : [200, 3500],
                            'nbins' : 50,
                            'label' : 'GenMuon_pt1 [GeV]',
                            'SafeCuts' : 'event.nMuons > 0',
}

variables['Muon2_genpt'] = {'name' : 'event.Muon_Genpt[1]',
                            'range' : [200, 3500],
                            'nbins' : 50,
                            'label' : 'GenMuon_pt2 [GeV]',
                            'SafeCuts' : 'event.nMuons > 1',
}


variables['Muon1_nGeomDets'] = {'name' : 'event.Muon_nGeomDets[0]',
                            'range' : [0, 50],
                            'nbins' : 50,
                            'label' : 'Muon_nCompGeomDets1',
                            'SafeCuts' : 'event.nMuons > 0',
}

variables['Muon2_nGeomDets'] = {'name' : 'event.Muon_nGeomDets[1]',
                            'range' : [0, 50],
                            'nbins' : 50,
                            'label' : 'Muon_nCompGeomDets2',
                            'SafeCuts' : 'event.nMuons > 1',
}


variables['Muon1_nHits'] = {'name' : 'event.Muon_nHits[0]',
                            'range' : [0, 170],
                            'nbins' : 50,
                            'label' : 'Muon_nCompHits1',
                            'SafeCuts' : 'event.nMuons > 0',
}

variables['Muon2_nHits'] = {'name' : 'event.Muon_nHits[1]',
                            'range' : [0, 170],
                            'nbins' : 50,
                            'label' : 'Muon_nCompHits2',
                            'SafeCuts' : 'event.nMuons > 1',
}

variables['Hit1_Muon1_GeomDet1_x'] = {'name' : 'event.Hit_x[0]',
                                      'range' : [-200, 200],
                                      'nbins' : 50,
                                      'label' : 'Hit_Muon1_GeomDet1_x [cm]',
                                      'SafeCuts' : '(event.nMuons > 0 and event.Muon_nGeomDets[0] > 0 and event.Muon_nHits[0] > 0)',
}

variables['Prop_Muon1_GeomDet1_x'] = {'name' : 'event.Prop_x[0]',
                                      'range' : [-200, 200],
                                      'nbins' : 50,
                                      'label' : 'Prop_Muon1_GeomDet1_x [cm]',
                                      'SafeCuts' : '(event.nMuons > 0 and event.Muon_nGeomDets[0] > 0)',
}




