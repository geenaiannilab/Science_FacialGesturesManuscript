
close all; 
clear all;

date = '171010';
subject ='Thor';
sess2run = 6;
statistic = 'variance';

% load bhv data
bhvDir = (['D:\Dropbox\Thor_FTExperiment\SUAinfo\' subject '_' date '/Data4Analysis']);
filePattern = fullfile(bhvDir, ['bhvspikes_sub' num2str(sess2run) '*']);
bhvFile = dir(filePattern);
TDTsessID = cell2mat(extractBetween(bhvFile.name,['sub' num2str(sess2run) '_'],['_']));
bhvFile = bhvFile.name;

% set up files to load 
STMAdir = ['D:\Dropbox\FaceExpressionAnalysis\STMA2021\' subject '_' date];
workDir = [STMAdir '\' subject '-' date '-' TDTsessID '\'];
f2load = dir([workDir 'cropped*']);
videoFileName = [f2load.folder '\' f2load.name];

correlationMatrix = correlatefROIs('subject', subject, 'videoFileName',videoFileName, 'pathToROIs', workDir,'statistic',statistic);


