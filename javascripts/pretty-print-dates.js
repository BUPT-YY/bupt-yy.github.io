$(function() {
    $('span.date').each( function() {
	var day = moment($(this).text(), 'YYYY-MM-DD hh:mm:ss a');
	$(this).text(moment(day).format('LLL'));
	$(this).attr('title', moment(day).fromNow());
    });
});
