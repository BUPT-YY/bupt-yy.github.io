function computeAnswer()
{
  var atsign = '@';

  var reg = /(J[a-z]*) (F[a-z]*)/;
  var matches = reg.exec( document.body.innerHTML );
  var lastName = matches[2].toLowerCase();

  var reg = /(math.[a-z]*.edu)/;
  var matches = reg.exec( document.body.innerHTML );
  var address = matches[1].toLowerCase();

  return lastName + atsign + address;
}

function sendLetter()
{
  var answer = computeAnswer();
  parent.location = 'm'+'a'+  'i'+  'l'+ 't' +  'o' + ':' + answer;
}

function phoneNumber()
{
  number = (6*6*23 + 1).toString();
  number = number + number.substring(1,2) + number.substring(0,1) + number.substring(2,3);
  number = (7 * 31 * 43 * parseInt(number)).toString();
  return '(' + number.substring(0,3) + ') ' + number.substring(3,6) + '-' + number.substring(6,10);
}

$(document).ready(function() {
  $('span.email').text( computeAnswer() );
  $('a.email').attr( 'href', 'm'+'a'+  'i'+  'l'+ 't' +  'o' + ':' + computeAnswer() );
  $('span.phone').text( phoneNumber() );
});
