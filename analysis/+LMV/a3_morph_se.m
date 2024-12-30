%% Morph align trials across recordings
%{
    NOTE:
    This is different from morphing each recording on their own. For a single recording, the target 
    is determined only by its own trials. Here, the target is determined by median intervals across 
    trials pooled from all recordings.

    5/17/24 run includes the following recordings
    
       recId                      folder                          name                                  path                            Subject     Block      Region         Subregion                  Coords                NativeCoords         LocalizationNote                                      LocalPathology                                 ECoGChan      StimNote           Task          NumTrials    Speech      FaceVideo                                                                                       KilosortPath                                                                                          Motion           SortTime         SortStatus       SurfaceCh_depth_nCh        Plots         Rating                                         Note                                           LG_ML         JC_KS   
    ___________    ____________________________________    __________________    ___________________________________________________    ________    ______    _________    _______________    _____________________________    ____________    __________________________    ________________________________________________________________________    _________    __________    ________________    _________    ______    _____________    ______________________________________________________________________________________________________________________________________________________________________________    _________________    ___________    ________________    ___________________    ______________    ______    ______________________________________________________________________________    __________    __________

    {'NP38_B5'}    {'C:\chang_lab\project_np\se\IFG'  }    {'NP38_B5_se.mat'}    {'C:\chang_lab\project_np\se\IFG\NP38_B5_se.mat'  }    {'NP38'}    {'B5'}    {'IFG'  }    {'opercularis'}    {0×0 char                   }        NaN         {'no brainlab, need help'}    {0×0 char                                                              }    {'NA'   }    {0×0 char}    {'LMV, SentGen'}       NaN       {'v3'}    {'captured' }    {'/data_store2/neuropixels/preproc/NP38_B5/kilosort/catgt_NP38_B5_g0/NP38_B5_g0_imec0/'                                                                                      }    {'stable'       }     14    1055    {'v2_mx'       }     {[NaN 6430 NaN]}      {0×0 char    }       4      {'site1_*.jpg'                                                               }    {0×0 char}    {0×0 char}
    {'NP38_B6'}    {'C:\chang_lab\project_np\se\IFG'  }    {'NP38_B6_se.mat'}    {'C:\chang_lab\project_np\se\IFG\NP38_B6_se.mat'  }    {'NP38'}    {'B6'}    {'IFG'  }    {'opercularis'}    {0×0 char                   }        NaN         {'no brainlab, need help'}    {0×0 char                                                              }    {'NA'   }    {0×0 char}    {'LMV'         }        97       {'v3'}    {'captured' }    {'/data_store2/neuropixels/preproc/NP38_B6/kilosort/catgt_NP38_B6_g0/NP38_B6_g0_imec0'                                                                                       }    {'stable'       }      0     607    {'v2_mx'       }     {[NaN 5830 NaN]}      {0×0 char    }       5      {'site2_*.jpg'                                                               }    {0×0 char}    {0×0 char}
    {'NP69_B1'}    {'C:\chang_lab\project_np\se\IFG'  }    {'NP69_B1_se.mat'}    {'C:\chang_lab\project_np\se\IFG\NP69_B1_se.mat'  }    {'NP69'}    {'B1'}    {'IFG'  }    {'opercularis'}    {[ -56.3000 33.2000 -7.2000]}        NaN         {0×0 char                }    {0×0 char                                                              }    {'NA'   }    {0×0 char}    {'LMV, semsr'  }       NaN       {'v3'}    {'captured' }    {'/data_store2/neuropixels/preproc/NP69_B1/kilosort/catgt_NP69_B1_g0/NP69_B1_g0_imec0, /data_store2/neuropixels/preproc/NP69_B1/kilosort/catgt_NP69_B1_g0/NP69_B1_g0_imec1'  }    {'stable;stable'}    173    1277    {'v2_mx; v2_mx'}     {[NaN 6925 NaN]}      {'yes'       }     NaN      {'Double probe'                                                              }    {0×0 char}    {0×0 char}
    {'NP69_B2'}    {'C:\chang_lab\project_np\se\IFG'  }    {'NP69_B2_se.mat'}    {'C:\chang_lab\project_np\se\IFG\NP69_B2_se.mat'  }    {'NP69'}    {'B2'}    {'IFG'  }    {'opercularis'}    {[ -56.3000 33.2000 -7.2000]}        NaN         {0×0 char                }    {0×0 char                                                              }    {'NA'   }    {0×0 char}    {'LMV, semsr'  }       NaN       {'v3'}    {'captured' }    {'/data_store2/neuropixels/preproc/NP69_B2/kilosort/catgt_NP69_B2_g0/NP69_B2_g0_imec0, /data_store2/neuropixels/preproc/NP69_B2/kilosort/catgt_NP69_B2_g0/NP69_B2_g0_imec1'  }    {'stable;stable'}     87    1029    {'v2_mx; v2_mx'}     {[NaN 7100 NaN]}      {'yes'       }     NaN      {'Double probe'                                                              }    {0×0 char}    {0×0 char}
    {'NP78_B1'}    {'C:\chang_lab\project_np\se\IFG'  }    {'NP78_B1_se.mat'}    {'C:\chang_lab\project_np\se\IFG\NP78_B1_se.mat'  }    {'NP78'}    {'B1'}    {'IFG'  }    {'opercularis'}    {[     -54.8000 34.7000 -16]}        NaN         {0×0 char                }    {'Scant infiltrating neoplastic astrocytes in subcortical white matter'}    {'NA'   }    {0×0 char}    {'LMV, semsr'  }       NaN       {'v3'}    {'captured' }    {'/data_store2/neuropixels/preproc/NP78_B1/kilosort/catgt_NP78_B1_g0/NP78_B1_g0_imec0'                                                                                       }    {'stable'       }    523    2156    {'v2_qg'       }     {[NaN 7550 NaN]}      {'yes'       }     NaN      {'semsr1 version'                                                            }    {0×0 char}    {0×0 char}
    {'NP35_B2'}    {'C:\chang_lab\project_np\se\STG'  }    {'NP35_B2_se.mat'}    {'C:\chang_lab\project_np\se\STG\NP35_B2_se.mat'  }    {'NP35'}    {'B2'}    {'mSTG' }    {0×0 char     }    {[ 73.2000 -4.3000 -13.4000]}        NaN         {0×0 char                }    {0×0 char                                                              }    {[  152]}    {0×0 char}    {'LMV, TIMIT*' }       NaN       {'v3'}    {'processed'}    {'/data_store2/neuropixels/preproc/NP35_B2/kilosort/mc_NP35_B2_g0/NP35_B2_g0_imec0'                                                                                          }    {'corrected'    }      0     Inf    {'v2_mx'       }     {[NaN 7300 NaN]}      {0×0 char    }       6      {'Used incorrect TIMIT from /ECOG_EXP/TIMIT/run_TIMIT_listen_short_intraop.m'}    {'Yes'   }    {0×0 char}
    {'NP43_B1'}    {'C:\chang_lab\project_np\se\STG'  }    {'NP43_B1_se.mat'}    {'C:\chang_lab\project_np\se\STG\NP43_B1_se.mat'  }    {'NP43'}    {'B1'}    {'pSTG' }    {0×0 char     }    {0×0 char                   }        NaN         {'need help'             }    {0×0 char                                                              }    {'NA'   }    {0×0 char}    {'LMV, TIMIT'  }       NaN       {'v3'}    {'captured' }    {'/data_store2/neuropixels/preproc/NP43_B1/kilosort/catgt_NP43_B1_g0/NP43_B1_g0_imec0/'                                                                                      }    {'ok'           }      0     Inf    {'v2_mx'       }     {[NaN 5730 NaN]}      {'yes (auto)'}       4      {'Probe degraded during recording. Added 'blike' to librispeech_lex.txt.'    }    {0×0 char}    {0×0 char}
    {'NP50_B3'}    {'C:\chang_lab\project_np\se\STG'  }    {'NP50_B3_se.mat'}    {'C:\chang_lab\project_np\se\STG\NP50_B3_se.mat'  }    {'NP50'}    {'B3'}    {'mSTG' }    {0×0 char     }    {[-61.1000 15.6000 -20.5000]}        NaN         {0×0 char                }    {0×0 char                                                              }    {'NA'   }    {0×0 char}    {'TIMIT, LMV'  }       NaN       {'v3'}    {'captured' }    {'/data_store2/neuropixels/preproc/NP50_B3/kilosort/mc_NP50_B3_g0/NP50_B3_g0_imec0'                                                                                          }    {'corrected'    }     80    1174    {'v2_mx'       }     {[NaN 6430 NaN]}      {0×0 char    }       2      {0×0 char                                                                    }    {0×0 char}    {0×0 char}
    {'NP53_B1'}    {'C:\chang_lab\project_np\se\STG'  }    {'NP53_B1_se.mat'}    {'C:\chang_lab\project_np\se\STG\NP53_B1_se.mat'  }    {'NP53'}    {'B1'}    {'aSTG' }    {0×0 char     }    {[-58.5000 19.8000 -25.3000]}        NaN         {0×0 char                }    {0×0 char                                                              }    {[   67]}    {0×0 char}    {'LMV, semsr'  }       NaN       {'v3'}    {'captured' }    {'/data_store2/neuropixels/preproc/NP53_B1/kilosort/catgt_NP53_B1_g0/NP53_B1_g0_imec0'                                                                                       }    {'stable'       }    102    1171    {'v2_mx'       }     {[NaN 6480 NaN]}      {'yes'       }       3      {0×0 char                                                                    }    {0×0 char}    {0×0 char}
    {'NP41_B1'}    {'C:\chang_lab\project_np\se\mPrCG'}    {'NP41_B1_se.mat'}    {'C:\chang_lab\project_np\se\mPrCG\NP41_B1_se.mat'}    {'NP41'}    {'B1'}    {'mPrCG'}    {0×0 char     }    {[ -49.5000 25.1000 30.3000]}        NaN         {0×0 char                }    {0×0 char                                                              }    {[86 87]}    {0×0 char}    {'LMV'         }        64       {'v3'}    {'NA'       }    {'/data_store2/neuropixels/preproc/NP41_B1/kilosort/NP41_B1_g0/NP41_B1_g0_imec0'                                                                                             }    {'stable'       }    200     730    {'v1.5_mx'     }     {[NaN 5870 NaN]}      {0×0 char    }       5      {0×0 char                                                                    }    {0×0 char}    {0×0 char}
    {'NP44_B2'}    {'C:\chang_lab\project_np\se\mPrCG'}    {'NP44_B2_se.mat'}    {'C:\chang_lab\project_np\se\mPrCG\NP44_B2_se.mat'}    {'NP44'}    {'B2'}    {'mPrCG'}    {0×0 char     }    {[      -54.9000 13 28.2000]}        NaN         {0×0 char                }    {0×0 char                                                              }    {'NA'   }    {0×0 char}    {'LMV'         }        72       {'v3'}    {'captured' }    {'/data_store2/neuropixels/preproc/NP44_B2/kilosort/mc_NP44_B2_g0/NP44_B2_g0_imec0/'                                                                                         }    {'corrected'    }     80    1077    {'v2_mx'       }     {[NaN 6720 NaN]}      {0×0 char    }       5      {0×0 char                                                                    }    {0×0 char}    {0×0 char}
    {'NP44_B3'}    {'C:\chang_lab\project_np\se\mPrCG'}    {'NP44_B3_se.mat'}    {'C:\chang_lab\project_np\se\mPrCG\NP44_B3_se.mat'}    {'NP44'}    {'B3'}    {'mPrCG'}    {0×0 char     }    {[ -53.1000 16.8000 30.6000]}        NaN         {0×0 char                }    {0×0 char                                                              }    {'NA'   }    {0×0 char}    {'LMV'         }        55       {'v3'}    {'captured' }    {'/data_store2/neuropixels/preproc/NP44_B3/kilosort/mc_NP44_B3_g0/NP44_B3_g0_imec0_ksdc-nb0/'                                                                                }    {'corrected'    }     82     780    {'v2_mx'       }     {[NaN 6810 NaN]}      {0×0 char    }       3      {0×0 char                                                                    }    {0×0 char}    {0×0 char}
    {'NP45_B1'}    {'C:\chang_lab\project_np\se\mPrCG'}    {'NP45_B1_se.mat'}    {'C:\chang_lab\project_np\se\mPrCG\NP45_B1_se.mat'}    {'NP45'}    {'B1'}    {'mPrCG'}    {0×0 char     }    {[       59 14.6000 30.1000]}        NaN         {0×0 char                }    {0×0 char                                                              }    {'NA'   }    {0×0 char}    {'LMV'         }        70       {'v3'}    {'captured' }    {'/data_store2/neuropixels/preproc/NP45_B1/kilosort/catgt_NP45_B1_g0/NP45_B1_g0_imec0/'                                                                                      }    {'stable'       }    118     685    {'v2_mx'       }     {[NaN 7440 NaN]}      {0×0 char    }       2      {0×0 char                                                                    }    {0×0 char}    {0×0 char}
    {'NP45_B2'}    {'C:\chang_lab\project_np\se\mPrCG'}    {'NP45_B2_se.mat'}    {'C:\chang_lab\project_np\se\mPrCG\NP45_B2_se.mat'}    {'NP45'}    {'B2'}    {'mPrCG'}    {0×0 char     }    {[       47.7000 10.3000 43]}        NaN         {0×0 char                }    {0×0 char                                                              }    {'NA'   }    {0×0 char}    {'LMV'         }        72       {'v3'}    {'captured' }    {'/data_store2/neuropixels/preproc/NP45_B2/kilosort/catgt_NP45_B2_g0/NP45_B2_g0_imec0/'                                                                                      }    {'corrected'    }     96     695    {'v2_mx'       }     {[NaN 6230 NaN]}      {0×0 char    }       4      {'Previously [NaN,7200,NaN]'                                                 }    {0×0 char}    {0×0 char}
    {'NP46_B1'}    {'C:\chang_lab\project_np\se\mPrCG'}    {'NP46_B1_se.mat'}    {'C:\chang_lab\project_np\se\mPrCG\NP46_B1_se.mat'}    {'NP46'}    {'B1'}    {'mPrCG'}    {0×0 char     }    {0×0 char                   }        NaN         {'no brainlab'           }    {0×0 char                                                              }    {'NA'   }    {0×0 char}    {'LMV'         }        69       {'v3'}    {'captured' }    {'/data_store2/neuropixels/preproc/NP46_B1/kilosort/catgt_NP46_B1_g0/NP46_B1_g0_imec0'                                                                                       }    {'stable'       }    260     762    {'v2_mx'       }     {[NaN 6180 NaN]}      {0×0 char    }       3      {'Previously [NaN,6840,NaN]'                                                 }    {0×0 char}    {0×0 char}
    {'NP47_B3'}    {'C:\chang_lab\project_np\se\mPrCG'}    {'NP47_B3_se.mat'}    {'C:\chang_lab\project_np\se\mPrCG\NP47_B3_se.mat'}    {'NP47'}    {'B3'}    {'mPrCG'}    {0×0 char     }    {[  61.2000 19.6000 24.2000]}        NaN         {0×0 char                }    {0×0 char                                                              }    {'NA'   }    {0×0 char}    {'LMV'         }        72       {'v3'}    {'captured' }    {'/data_store2/neuropixels/preproc/NP47_B3/kilosort/catgt_NP47_B3_g0/NP47_B3_g0_imec0/'                                                                                      }    {'stable'       }      0     650    {'v1.5_mx'     }     {[NaN 6550 NaN]}      {0×0 char    }       2      {'Originally located as vPrCG.'                                              }    {0×0 char}    {0×0 char}
    {'NP56_B1'}    {'C:\chang_lab\project_np\se\mPrCG'}    {'NP56_B1_se.mat'}    {'C:\chang_lab\project_np\se\mPrCG\NP56_B1_se.mat'}    {'NP56'}    {'B1'}    {'mPrCG'}    {0×0 char     }    {[       59 14.6000 30.1000]}        NaN         {0×0 char                }    {0×0 char                                                              }    {'NA'   }    {0×0 char}    {'LMV'         }       NaN       {'v3'}    {'captured' }    {'/data_store2/neuropixels/preproc/NP56_B1/kilosort/catgt_NP56_B1_g0/NP56_B1_g0_imec0/; /data_store2/neuropixels/preproc/NP56_B1/kilosort/catgt_NP56_B1_g0/NP56_B1_g0_imec1/'}    {'ok;ok'        }    160    1576    {'v2_mx; v2_mx'}     {[NaN 6290 NaN]}      {0×0 char    }       7      {0×0 char                                                                    }    {0×0 char}    {0×0 char}
    {'NP52_B1'}    {'C:\chang_lab\project_np\se\vPrCG'}    {'NP52_B1_se.mat'}    {'C:\chang_lab\project_np\se\vPrCG\NP52_B1_se.mat'}    {'NP52'}    {'B1'}    {'vPrCG'}    {0×0 char     }    {[  -61.7000 22.3000 0.8000]}        NaN         {0×0 char                }    {0×0 char                                                              }    {'NA'   }    {0×0 char}    {'LMV'         }        72       {'v3'}    {'captured' }    {'/data_store2/neuropixels/preproc/NP52_B1/kilosort/catgt_NP52_B1_g0/NP52_B1_g0_imec0/'                                                                                      }    {'ok'           }    128     785    {'v2_mx'       }     {[NaN 6850 NaN]}      {0×0 char    }       4      {0×0 char                                                                    }    {0×0 char}    {0×0 char}
    {'NP52_B2'}    {'C:\chang_lab\project_np\se\vPrCG'}    {'NP52_B2_se.mat'}    {'C:\chang_lab\project_np\se\vPrCG\NP52_B2_se.mat'}    {'NP52'}    {'B2'}    {'vPrCG'}    {0×0 char     }    {[  -61.3000 22.9000 3.3000]}        NaN         {0×0 char                }    {0×0 char                                                              }    {'NA'   }    {0×0 char}    {'LMV'         }        72       {'v3'}    {'captured' }    {'/data_store2/neuropixels/preproc/NP52_B2/kilosort/catgt_NP52_B2_g0/NP52_B2_g0_imec0'                                                                                       }    {'ok'           }     36     810    {'v2_mx'       }     {[NaN 7000 NaN]}      {0×0 char    }       3      {'Surface not very clear.'                                                   }    {0×0 char}    {0×0 char}
    {'NP54_B1'}    {'C:\chang_lab\project_np\se\vPrCG'}    {'NP54_B1_se.mat'}    {'C:\chang_lab\project_np\se\vPrCG\NP54_B1_se.mat'}    {'NP54'}    {'B1'}    {'vPrCG'}    {0×0 char     }    {[       -61 20.2000 2.4000]}        NaN         {0×0 char                }    {0×0 char                                                              }    {'NA'   }    {0×0 char}    {'LMV'         }        72       {'v3'}    {'captured' }    {'/data_store2/neuropixels/preproc/NP54_B1/kilosort/catgt_NP54_B1_g0/NP54_B1_g0_imec0'                                                                                       }    {'ok'           }    250     850    {'v2_mx'       }     {[NaN 6670 NaN]}      {0×0 char    }       3      {0×0 char                                                                    }    {0×0 char}    {0×0 char}

%}

anaDir = LMV.Data.GetAnalysisDir('data');
srcTb = LMV.Data.FindSource([]);

%% Morph align repetitions of the same sentences across all recordings (M1 for short)

m1Dir = fullfile(anaDir, 'se_m1');
if ~exist(m1Dir, 'dir')
    mkdir(m1Dir);
end

m1Paths = fullfile(m1Dir, strrep(srcTb.name, '_se.mat', '_se_m1.mat'));

if all(cellfun(@(x) exist(x,'file'), m1Paths))
    % Load cached results
    seM1Array = NP.SE.LoadSession(m1Paths);
    
else
    % Load recordings
    seArray = NP.SE.LoadSession(srcTb.path, 'UserFunc', @(x) x.RemoveTable('LFP', 'niTime'));
    
    % Enrich se
    ops = NP.Param.Enrich;
    ops.isMergeMeta = true;
    ops.isSpkRate = true;
    ops.isMel = true;
    ops.isPitch = true;
    ops.isArtic = true;
    for i = 1 : numel(seArray)
        seArray(i) = LMV.SE.Transform(seArray(i), 'enrich', ops);
    end
    
    % Concatenate taskTime and taskValue tables across recordings
    seTaskArray = seArray.Duplicate({'taskTime', 'taskValue'});
    seTaskCat = NP.SE.CatEpochs(seTaskArray);
    
    % Finding morph times
    NP.SE.SetMorphTimes(seTaskCat);
    
    % Split concatenated taskTime and taskValue tables back to individual recordings
    seTaskArray = seTaskCat.Split([seTaskArray.numEpochs]');
    
    for i = numel(seArray) : -1 : 1
        % Set updated tables and reference times back to the original se
        seArray(i).SetTable('taskTime', seTaskArray(i).GetTable('taskTime'));
        seArray(i).SetTable('taskValue', seTaskArray(i).GetTable('taskValue'));
        seArray(i).SetReferenceTime(seTaskArray(i).GetReferenceTime);
        
        % Run morphing on the recording
        seM1Array(i,1) = NP.SE.MorphSession(seArray(i));
        
        % Make sure all spike times remain in cell array
        seM1Array(i,1).Column2Cell('spikeTime');
    end
    
    % Make sure all prod events remain in cell arrays for later concatenation
    LMV.SE.StandardizeDataType(seM1Array);
    
    % Cache morphed se
    for i = 1 : numel(seM1Array)
        se = seM1Array(i);
        save(m1Paths{i}, 'se', '-v7.3');
    end
end

%% Further morph align trials by LMV task events (M2 for short)

m2Dir = fullfile(anaDir, 'se_m2');
m2Paths = fullfile(m2Dir, strrep(srcTb.name, '_se.mat', '_se_m2.mat'));

if all(cellfun(@(x) exist(x,'file'), m2Paths))
    % Load cached results
    seM2Array = NP.SE.LoadSession(m2Paths);
    
else
    % Load recordings
    if ~exist('seM1Array', 'var')
        seM1Array = NP.SE.LoadSession(m1Paths);
    end
    
    % Concatenate taskTime and taskValue tables across recordings
    seM1TaskArray = seM1Array.Duplicate({'taskTime', 'taskValue'});
    seM1Cat = NP.SE.CatEpochs(seM1TaskArray);
    
    % Finding morph times
    NP.SE.SetMorphTimes(seM1Cat, 'events', {'cue1On', 'cue1Off', 'stimOn', 'stimOff', 'cue3On', 'prodMatchOn', 'prodMatchOff'});
    
    % Split concatenated taskTime and taskValue tables back to individual recording's
    seM1TaskArray = seM1Cat.Split([seM1TaskArray.numEpochs]');
    
    for i = numel(seM1Array) : -1 : 1
        % Set updated tables and reference times back to the original se
        seM1Array(i).SetTable('taskTime', seM1TaskArray(i).GetTable('taskTime'));
        seM1Array(i).SetTable('taskValue', seM1TaskArray(i).GetTable('taskValue'));
        seM1Array(i).SetReferenceTime(seM1TaskArray(i).GetReferenceTime);
        
        % Run morphing on the recording
        seM2Array(i,1) = NP.SE.MorphSession(seM1Array(i));
    end
    
    % Make sure all prod events remain in cell arrays for later concatenation
    LMV.SE.StandardizeDataType(seM2Array);
    
    % Cache morphed se
    if ~exist(m2Dir, 'dir')
        mkdir(m2Dir);
    end
    for i = 1 : numel(seM2Array)
        se = seM2Array(i);
        save(m2Paths{i}, 'se', '-v7.3');
    end
end

%% Prepare data for plotting

% Load original
if ~exist('seArray', 'var')
    seArray = NP.SE.LoadSession(srcTb.path, ...
        'UserFunc', @(x) x.RemoveTable('LFP', 'niTime', 'spikeTime'));
    ops = NP.Param.Enrich;
    ops.isMel = true;
    for i = 1 : numel(seArray)
        seArray(i) = LMV.SE.Transform(seArray(i), 'enrich', ops);
    end
end

% Make lite se with only the required data
seAll = [seArray seM1Array seM2Array];
seAll = seAll.Duplicate({'taskTime', 'taskValue', 'mel'});

% Concatenate recordings
for i = size(seAll,2) : -1 : 1
    seCats(i) = Merge(seAll(:,i));
    seCats(i).userData = seCats(i).userData(1);
end

% Remove bad trials
isBad = LMV.SE.IsBadTrials(seCats(2), -Inf, 0.5);
arrayfun(@(x) x.RemoveEpochs(isBad), seCats);

senTbs = cell(size(seCats));
for i = 1 : numel(senTbs)
    % Split and group trials by sentence
    senTb = NP.TaskBaseClass.SplitBySentence(seCats(i));
    
    % Align to prodMatchOn or prodOn
    for j = 1 : height(senTb)
        tt = senTb.se(j).GetTable('taskTime');
        if i == 1
            t0 = cellfun(@(x) x(1), tt.prodOn); % for original se
        else
            t0 = median(tt.prodMatchOn); % for the morphed
        end
        senTb.se(j).AlignTime(t0);
    end
    
    senTbs{i} = senTb;
end

%% Plot alignments of each sentence

figDir = fullfile(anaDir, 'alignment');
if ~exist(figDir, "dir")
    mkdir(figDir);
end

senInd = 1 : height(senTbs{1});
tWin = [-0.5 2.5];

for i = 1 : numel(senInd)
    f = MPlot.Figure(10); clf
    f.Position = [52 400 1800 500];
    
    rowDist = [1 3];
    colDist = ones(size(senTbs));
    
    tl = tiledlayout(sum(rowDist), sum(colDist));
    tl.Padding = "compact";
    
    for j = 1 : numel(senTbs)
        % Get data
        se = senTbs{j}.se(senInd(i));
        tv = se.GetTable('taskValue');
        
        % Plot phonetic labels of the first trial
        ntArgs = MPlot.FindTileInd(rowDist, colDist, 1, j);
        ax = nexttile(ntArgs{:});
        NP.TaskBaseClass.PlotTGEHier(ax, se, 1, [], 'prod');
        ax.XLim = tWin;
        ax.XAxis.Visible = 'off';
        
        % Plot phoneme patches of all trials
        ntArgs = MPlot.FindTileInd(rowDist, colDist, 2, j);
        ax = nexttile(ntArgs{:});
        NP.TaskBaseClass.PlotTGE(ax, se, [], [-Inf Inf], 'phone', 'prod', 'Style', 'patch', 'FontSize', 0);
        
        wRib = [-0.5 -0.45];
        cFunc = @(x) brewermap(x, 'Set3');
        [xSeg, ySeg] = MPlot.GroupRibbon(wRib, tv.recId, cFunc, 'Style', 'patch', 'PlotArgs', {'Parent', ax});
        text(xSeg-.1, sort(ySeg), unique(tv.recId, 'stable'), 'Hori', 'right', 'Interpreter', 'none');
        
        ax.XLim = tWin;
        ax.YAxis.Visible = 'off';
        ax.Title.String = sprintf("N = %i trials", se.numEpochs);
    end
    
    senIdx = find(strcmp(tv.stimId(1), LMV.Param.stimIdList), 1);
    
    exportgraphics(f, fullfile(figDir, sprintf("morphing_sent%i_%s.png", senIdx, tv.stimId(1))));
end

%% Plot alignment of all sentences

tWin = [-0.5 2.5];

f = MPlot.Figure(11); clf
f.Position = [52 589 1800 600];

rowDist = 1;
colDist = ones(size(senTbs));

tl = tiledlayout(sum(rowDist), sum(colDist));
tl.Padding = "compact";

for j = 1 : numel(senTbs)
    % Exclude less repeated sentences from a few subjects with early task versions
    senTb = senTbs{j};
    isFew = senTb.numTrial < 20;
    senTb(isFew,:) = [];
    
    % Combine se objects across sentences
    se = Merge(senTb.se);
    tv = se.GetTable('taskValue');
    stimText = cellstr(senTb.stimText);
    stimText = cellfun(@(x) [x(1:10) '...'], stimText, 'Uni', false);
    
    % Plot phoneme patches of all trials
    ntArgs = MPlot.FindTileInd(rowDist, colDist, 1, j);
    ax = nexttile(ntArgs{:});
    NP.TaskBaseClass.PlotTGE(ax, se, [], [-Inf Inf], 'phone', 'prod', 'Style', 'patch', 'FontSize', 0);
    
    wRib = tWin(1) + [.05 0.1];
    cFunc = @(x) brewermap(x, 'Set3');
    [xSeg, ySeg] = MPlot.GroupRibbon(wRib, tv.recId, cFunc, 'Style', 'patch', 'PlotArgs', {'Parent', ax});
    
    wRib = tWin(1) + [0 .05];
    cFunc = @(x) brewermap(x, 'Spectral');
    [xSeg, ySeg] = MPlot.GroupRibbon(wRib, tv.stimId, cFunc, 'Style', 'patch', 'PlotArgs', {'Parent', ax});
    text(xSeg-.1, sort(ySeg), stimText, 'Hori', 'right', 'Interpreter', 'none');
    
    ax.XLim = tWin;
    ax.YAxis.Visible = 'off';
    ax.Title.String = sprintf("N = %i trials", se.numEpochs);
end

exportgraphics(f, fullfile(figDir, "morphing_all_sent.png"));

%% Plot alignment of all task events (sentence only)

tWin = [-.5 7.5];

f = MPlot.Figure(12); clf
f.Position = [52 589 1800 600];

rowDist = 1;
colDist = ones(size(senTbs));

tl = tiledlayout(sum(rowDist), sum(colDist));
tl.Padding = "compact";

for j = 1 : numel(senTbs)
    % Exclude less repeated sentences from a few subjects with early task versions
    senTb = senTbs{j};
    isFew = senTb.numTrial < 20;
    senTb(isFew,:) = [];
    
    % Combine se objects across sentences
    se = Merge(senTb.se);
    se.AlignTime('cue1On', 'taskTime');
    tv = se.GetTable('taskValue');
    stimText = cellstr(senTb.stimText);
    stimText = cellfun(@(x) [x(1:10) '...'], stimText, 'Uni', false);
    
    % Plot task events of all trials
    ntArgs = MPlot.FindTileInd(rowDist, colDist, 1, j);
    ax = nexttile(ntArgs{:});
    colNames = {'cue1', 'cue3', 'prod', 'stim'};
    for k = 1 : se.numEpochs
        NP.TaskBaseClass.PlotEventWindows(ax, se, k, tWin, colNames, 'Style', 'block', 'Alpha', .5);
    end
    
    % Plot ribbon labels
    wRib = tWin(1) + [.05 0.1];
    cFunc = @(x) brewermap(x, 'Set3');
    [xSeg, ySeg] = MPlot.GroupRibbon(wRib, tv.recId, cFunc, 'Style', 'patch', 'PlotArgs', {'Parent', ax});
    
    wRib = tWin(1) + [0 .05];
    cFunc = @(x) brewermap(x, 'Spectral');
    [xSeg, ySeg] = MPlot.GroupRibbon(wRib, tv.stimId, cFunc, 'Style', 'patch', 'PlotArgs', {'Parent', ax});
    text(xSeg-.1, sort(ySeg), stimText, 'Hori', 'right', 'Interpreter', 'none');
    
    ax.XLim = tWin;
    ax.YLim = [0.5 se.numEpochs+.5];
    ax.YAxis.Visible = 'off';
    ax.Title.String = sprintf("N = %i trials", se.numEpochs);
end

exportgraphics(f, fullfile(figDir, "morphing_task.png"));

%% Plot alignment of all task events (including individual phonemes)

tWin = [-.5 7.5];

f = MPlot.Figure(13); clf
f.Position = [52 589 1800 600];

rowDist = 1;
colDist = ones(size(senTbs));

tl = tiledlayout(sum(rowDist), sum(colDist));
tl.Padding = "compact";

for j = 1 : numel(senTbs)
    % Exclude less repeated sentences from a few subjects with early task versions
    senTb = senTbs{j};
    isFew = senTb.numTrial < 20;
    senTb(isFew,:) = [];
    
    % Combine se objects across sentences
    se = Merge(senTb.se);
    se.AlignTime('cue1On', 'taskTime');
    tv = se.GetTable('taskValue');
    stimText = cellstr(senTb.stimText);
    stimText = cellfun(@(x) [x(1:10) '...'], stimText, 'Uni', false);
    
    % Plot phoneme patches of all trials
    ntArgs = MPlot.FindTileInd(rowDist, colDist, 1, j);
    ax = nexttile(ntArgs{:});
    NP.TaskBaseClass.PlotTGE(ax, se, [], [-Inf Inf], 'phone', {'prod', 'stim'}, 'Style', 'patch', 'FontSize', 0);
    
    % Plot cues
    colNames = {'cue1', 'cue3'};
    for k = 1 : se.numEpochs
        NP.TaskBaseClass.PlotEventWindows(ax, se, k, tWin, colNames, 'Style', 'block', 'Alpha', .5);
    end
    
    % Plot ribbon labels
    wRib = tWin(1) + [.05 0.1];
    cFunc = @(x) brewermap(x, 'Set3');
    [xSeg, ySeg] = MPlot.GroupRibbon(wRib, tv.recId, cFunc, 'Style', 'patch', 'PlotArgs', {'Parent', ax});
    
    wRib = tWin(1) + [0 .05];
    cFunc = @(x) brewermap(x, 'Spectral');
    [xSeg, ySeg] = MPlot.GroupRibbon(wRib, tv.stimId, cFunc, 'Style', 'patch', 'PlotArgs', {'Parent', ax});
    text(xSeg-.1, sort(ySeg), stimText, 'Hori', 'right', 'Interpreter', 'none');
    
    ax.XLim = tWin;
    ax.YLim = [0.5 se.numEpochs+.5];
    ax.YAxis.Visible = 'off';
    ax.Title.String = sprintf("N = %i trials", se.numEpochs);
end

exportgraphics(f, fullfile(figDir, "morphing_task_detailed.png"));
