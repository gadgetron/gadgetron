function configure_accumulator(XMLconfig)
   global accumulator_buffer
   global accumulator_center_line 

   matrix_size = [str2num(XMLGetXPath(XMLconfig, '//ismrmrdHeader/encoding/reconSpace/matrixSize/x')), ...
			 str2num(XMLGetXPath(XMLconfig, '//ismrmrdHeader/encoding/reconSpace/matrixSize/y')), ...
			 str2num(XMLGetXPath(XMLconfig, '//ismrmrdHeader/encoding/reconSpace/matrixSize/z')), ...
                         str2num(XMLGetXPath(XMLconfig, '//ismrmrdHeader/acquisitionSystemInformation/receiverChannels'))];


   accumulator_center_line = str2num(XMLGetXPath(XMLconfig, '//ismrmrdHeader/encoding/encodingLimits/kspace_encoding_step_1/center'));

   fprintf('Accumulator: Reconstructing on matrix [%d, %d, %d]\n', matrix_size(1), matrix_size(2), matrix_size(3));


   accumulator_buffer = single(zeros(matrix_size));
   
end
